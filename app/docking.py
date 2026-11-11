# ══ CELL 10 · app/docking.py ══
"""
docking.py — High-level docking pipeline (SMILES → engine → result dict).

    SMILES → RDKit 3D mol → PDBQT (Meeko) + SDF → chosen engine → pose + energy

Also owns the runtime registry of targets (receptor paths + search boxes).
"""
import logging
from pathlib import Path
from typing import Any, Dict, Optional

from rdkit import Chem

from .chem import compute_qed, mol_properties, mol_to_pdbqt_string, mol_to_sdf_block, smiles_to_mol3d
from .lipinski import lipinski_filter, rejection_message
from .synth import synth_filter, synth_rejection_message
from .config import MAX_BOX_EDGE, MIN_BOX_EDGE, RECEPTORS_DIR
from .engines import get_engine
from .engines.base import EngineError
from .engines.poses import pose_to_pdb
from .zero_gpu import engine_uses_gpu, run_engine_on_gpu

logger = logging.getLogger(__name__)

# Filled at startup by main.py once receptors are prepared
RECEPTOR_PATHS: Dict[str, Optional[Path]] = {}

# Search boxes are registered when the user imports or uploads a PDB.
BOX_PARAMS = {}


def clamp_box(center, size):
    c = tuple(float(v) for v in center)
    s = tuple(max(MIN_BOX_EDGE, min(MAX_BOX_EDGE, float(v))) for v in size)
    return c, s


def register_custom_target(target_id: str, pdbqt_path: Path, center, size):
    c, s = clamp_box(center, size)
    RECEPTOR_PATHS[target_id] = pdbqt_path
    BOX_PARAMS[target_id] = {"center": c, "size": s}


def set_box(target_id: str, center, size):
    c, s = clamp_box(center, size)
    BOX_PARAMS[target_id] = {"center": c, "size": s}
    return c, s


def receptor_pdb_for(pdbqt_path: Optional[Path]) -> Optional[Path]:
    """Protonated PDB that sits next to the receptor PDBQT (<name>_H.pdb)."""
    if not pdbqt_path:
        return None
    cand = Path(pdbqt_path).with_name(Path(pdbqt_path).stem + "_H.pdb")
    return cand if cand.exists() else None


def dock_smiles(
    target_id: str,
    smiles: str,
    exhaustiveness: int = 8,
    engine_name: str = "vina",
    lipinski_mode: str = "standard",
    synth_mode: str = "standard",
    box_override: Optional[Dict[str, Any]] = None,
    n_poses: int = 10,
    diffdock_api_key: Optional[str] = None,
) -> Dict[str, Any]:
    """Full pipeline. Returns a dict (see keys below). Raises on any failure — no mock results."""
    if target_id not in BOX_PARAMS:
        raise ValueError(f"Unknown target: {target_id!r}. Choose from {list(BOX_PARAMS)}")

    smiles = (smiles or "").strip()
    if not smiles:
        raise ValueError("SMILES is empty — load or draw a molecule first.")

    # Lipinski gate: cheap 2-D check before 3-D embedding, Meeko or docking.
    lipinski_mode = (lipinski_mode or "standard").strip().lower()
    base_mol = Chem.MolFromSmiles(smiles)
    if base_mol is None:
        raise ValueError(f"Invalid SMILES: {smiles!r}")
    lipinski1 = lipinski_filter(base_mol, mode=lipinski_mode)
    if not lipinski1["passed"]:
        raise ValueError(rejection_message(lipinski1))

    # Calculate SA for every ligand; "off" disables only the gate, not the score.
    synth_mode = (synth_mode or "standard").strip().lower()
    synth = synth_filter(base_mol, mode=synth_mode)
    if not synth["passed"]:
        raise ValueError(synth_rejection_message(synth))

    import os as _os, random as _random, sys as _sys
    _seed = int(_os.environ["MEDCHEM_SEED"]) if _os.environ.get("MEDCHEM_SEED") else _random.randint(1, 2**31 - 1)
    for _m in ("app.chem", "app.engines.vina_engine", "app.engines.gnina_engine",
               "app.engines.karmadock_engine", "app.engines.diffdock_engine"):
        _mod = _sys.modules.get(_m)
        if _mod is not None and hasattr(_mod, "RANDOM_SEED"):
            _mod.RANDOM_SEED = _seed
    mol = smiles_to_mol3d(smiles)                    # raises ValueError for bad SMILES
    qed_score = compute_qed(mol)
    props = mol_properties(mol)
    props["lipinski1"] = lipinski1
    ligand_pdbqt = mol_to_pdbqt_string(mol)
    ligand_sdf = mol_to_sdf_block(mol)

    engine = get_engine(engine_name)
    box = box_override or BOX_PARAMS[target_id]
    box = {
        "center": tuple(float(v) for v in box["center"]),
        "size": tuple(float(v) for v in box["size"]),
    }
    box_adjustment = None
    if engine.name == "vina":
        # Vina's grid memory grows with box volume. Keep its effective box within
        # the documented 30 Å per-edge safe search range for laptop use.
        requested_size = box["size"]
        effective_size = tuple(min(30.0, edge) for edge in requested_size)
        if effective_size != requested_size:
            box_adjustment = {
                "requested_size": list(requested_size),
                "used_size": list(effective_size),
                "reason": "Vina grid memory limit (30 Å maximum per edge)",
            }
            logger.warning(
                "Vina search box reduced from %s Å to %s Å per edge to avoid grid-memory allocation failure",
                requested_size, effective_size,
            )
            box["size"] = effective_size
    receptor_path = RECEPTOR_PATHS.get(target_id)

    logger.info("Docking %s → %s with %s (exh=%d)", smiles[:50], target_id, engine.display_name, exhaustiveness)
    def run_engine(search_size):
        engine_kwargs = {
            "receptor_path": receptor_path,
            "ligand_pdbqt": ligand_pdbqt,
            "center": box["center"],
            "size": search_size,
            "exhaustiveness": exhaustiveness,
            "n_poses": max(1, min(20, int(n_poses or 10))),
            "smiles": smiles,
            "ligand_sdf": ligand_sdf,
            "receptor_pdb": receptor_pdb_for(receptor_path),
        }
        if engine.name == "diffdock":
            # This is request-scoped and is never written to process env or disk.
            engine_kwargs["diffdock_api_key"] = diffdock_api_key
        if engine_uses_gpu(engine.name):
            try:
                return run_engine_on_gpu(engine.name, engine_kwargs)
            except Exception as exc:
                if engine.name not in {"gnina", "karmadock"}:
                    raise
                # Preserve availability if the hosted GPU lease aborts or is unavailable.
                logger.warning("%s GPU execution failed; retrying on CPU: %s", engine.display_name, exc)
                engine_kwargs["force_cpu"] = True
                engine_kwargs["gpu_fallback_reason"] = f"GPU execution failed: {str(exc)[:250]}"
                return engine.dock(**engine_kwargs)
        return engine.dock(**engine_kwargs)

    try:
        result = run_engine(box["size"])
    except EngineError as exc:
        if engine.name != "vina" or "insufficient memory" not in str(exc).lower():
            raise
        last_error = exc
        original_size = box["size"]
        attempted_sizes = {original_size}
        # Grid allocation can still fail on low-memory laptops at 30 Å. Retry
        # progressively smaller boxes around the same binding-site centre.
        for factor in (0.8, 0.6, 0.45):
            retry_size = tuple(round(max(min(12.0, edge), edge * factor), 1) for edge in original_size)
            if retry_size in attempted_sizes:
                continue
            attempted_sizes.add(retry_size)
            try:
                result = run_engine(retry_size)
                if box_adjustment is None:
                    box_adjustment = {"requested_size": list(original_size)}
                box_adjustment.update({
                    "used_size": list(retry_size),
                    "reason": "Vina grid memory retry (smaller search box)",
                })
                box["size"] = retry_size
                logger.warning("Vina grid memory retry succeeded with box %s Å", retry_size)
                break
            except EngineError as retry_exc:
                if "insufficient memory" not in str(retry_exc).lower():
                    raise
                last_error = retry_exc
        else:
            raise EngineError(
                "Vina could not allocate its search grid even after reducing the box to a laptop-safe size. "
                "Reduce the box further or close other memory-heavy applications."
            ) from last_error
    logger.info("Done (%s): %.2f kcal/mol, QED %.3f", engine.display_name, result.energy, qed_score)

    return {
        "energy": result.energy,
        "qed": qed_score,
        "props": props,
        "synth": synth,
        "docked_pdbqt": result.docked_pdbqt,
        "docked_pdb": result.docked_pdb or pose_to_pdb(result.docked_pdbqt),
        "ligand_pdbqt": ligand_pdbqt,
        "smiles": smiles,
        "target_id": target_id,
        "engine": engine.name,
        "engine_metadata": {
            **result.metadata,
            "seed": _seed,
            **({"box_adjustment": box_adjustment} if box_adjustment else {}),
        },
        "is_mock": False,
        "box": {"center": list(box["center"]), "size": list(box["size"])},
    }


# ── CLI self-test:  python -m app.docking <loaded-target-id> vina ──────────────
if __name__ == "__main__":
    import sys
    logging.basicConfig(level=logging.INFO, format="%(levelname)s: %(message)s")
    from .receptor_prep import RECEPTORS
    for tid, info in RECEPTORS.items():
        p = RECEPTORS_DIR / f"{tid}_{info['pdb_id']}.pdbqt"
        RECEPTOR_PATHS[tid] = p if p.exists() else None
    if len(sys.argv) < 2:
        raise SystemExit("Import or upload a PDB in MobiDock first, then pass its target ID.")
    tid = sys.argv[1]
    eng = sys.argv[2] if len(sys.argv) > 2 else "vina"
    r = dock_smiles(tid, "CC(=O)Oc1ccccc1C(=O)O", engine_name=eng)
    print(f"energy={r['energy']}  qed={r['qed']}  engine={r['engine']}")
