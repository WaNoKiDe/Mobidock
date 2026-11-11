# ══ CELL 13 · app/engines/rescoring.py ══
"""
rescoring.py — Turn any 3D pose into a Vina-scale binding energy (kcal/mol).

ML engines (DiffDock-L, KarmaDock) output confidence/scores, not kcal/mol.
To keep the leaderboard fair we locally relax the predicted pose with the
AutoDock Vina scoring function (`vina --local_only`) and use THAT energy.
The engine's own score is kept in the result metadata.
"""
import logging
import re
import subprocess
import tempfile
from pathlib import Path
from typing import Optional, Tuple

from rdkit import Chem

from ..chem import mol_to_pdbqt_string
from ..config import find_vina, RANDOM_SEED
from .base import EngineError, EngineNotReady
from .poses import pose_centroid

logger = logging.getLogger(__name__)


def vina_rescoring_available() -> bool:
    """True when either the Vina CLI or its Python API can score poses."""
    if find_vina() is not None:
        return True
    try:
        from vina import Vina  # noqa: F401
        return True
    except Exception:
        return False


def _parse_vina_table(stdout: str) -> Optional[float]:
    for line in stdout.splitlines():
        m = re.match(r"\s+1\s+(-?\d+\.\d+)", line)
        if m:
            return float(m.group(1))
    m = re.search(r"Estimated Free Energy of Binding\s*:\s*(-?\d+\.\d+)", stdout)
    return float(m.group(1)) if m else None


def ensure_hs(mol: Chem.Mol) -> Chem.Mol:
    """ML engines often return heavy-atom-only poses; Meeko needs explicit hydrogens."""
    if any(a.GetAtomicNum() == 1 for a in mol.GetAtoms()):
        return mol
    try:
        Chem.SanitizeMol(mol)
    except Exception:
        pass
    return Chem.AddHs(mol, addCoords=True)


def centroid_in_box(mol: Chem.Mol, center, size, margin: float = 4.0) -> bool:
    conf = mol.GetConformer()
    atoms = [i for i, atom in enumerate(mol.GetAtoms()) if atom.GetAtomicNum() > 1]
    if not atoms:
        return False
    c = [sum(conf.GetAtomPosition(i)[k] for i in atoms) / len(atoms) for k in range(3)]
    return all(abs(c[k] - center[k]) <= size[k] / 2 + margin for k in range(3))


def vina_relax_and_score(
    receptor_pdbqt: Path,
    pose_mol: Chem.Mol,
    box_center: Tuple[float, float, float],
    box_size: Tuple[float, float, float],
    timeout: int = 180,
) -> Tuple[float, str]:
    """
    Locally optimise `pose_mol` (3D, with Hs) inside the receptor and score it.
    Returns (energy_kcal_per_mol, relaxed_pose_pdbqt).
    The search box is centred on the pose and enlarged so it always fits.
    """
    if receptor_pdbqt is None or not Path(receptor_pdbqt).exists():
        raise EngineNotReady("receptor PDBQT not prepared")

    pose_mol = ensure_hs(pose_mol)
    lig_pdbqt = mol_to_pdbqt_string(pose_mol)
    centroid = pose_centroid(lig_pdbqt) or box_center

    # box: centred on pose, at least as large as the game's box and the pose + margin
    conf = pose_mol.GetConformer()
    xs = [conf.GetAtomPosition(i).x for i in range(pose_mol.GetNumAtoms())]
    ys = [conf.GetAtomPosition(i).y for i in range(pose_mol.GetNumAtoms())]
    zs = [conf.GetAtomPosition(i).z for i in range(pose_mol.GetNumAtoms())]
    ext = [max(xs) - min(xs), max(ys) - min(ys), max(zs) - min(zs)]
    size = [max(box_size[i], ext[i] + 14.0) for i in range(3)]

    vina = find_vina()
    with tempfile.TemporaryDirectory() as td:
        lp, op = Path(td, "pose.pdbqt"), Path(td, "relaxed.pdbqt")
        lp.write_text(lig_pdbqt, encoding="utf-8")
        cli_error = None
        if vina is not None:
            cmd = [
                str(vina), "--receptor", str(receptor_pdbqt), "--ligand", str(lp), "--out", str(op),
                "--local_only",
                "--center_x", f"{centroid[0]:.3f}", "--center_y", f"{centroid[1]:.3f}", "--center_z", f"{centroid[2]:.3f}",
                "--size_x", f"{size[0]:.1f}", "--size_y", f"{size[1]:.1f}", "--size_z", f"{size[2]:.1f}",
            ]
            try:
                res = subprocess.run(cmd, capture_output=True, text=True, timeout=timeout)
                if res.returncode != 0:
                    raise EngineError((res.stderr or res.stdout)[-300:].strip() or f"Vina exited {res.returncode}")
                energy = _parse_vina_table(res.stdout)
                if energy is None:
                    raise EngineError("Could not read the Vina rescoring energy.")
                relaxed = op.read_text(encoding="utf-8") if op.exists() else lig_pdbqt
                return energy, relaxed
            except Exception as exc:
                cli_error = exc
                logger.warning("Vina CLI pose rescoring failed; trying the installed Python API: %s", exc)

        try:
            from vina import Vina
        except Exception as exc:
            if vina is None:
                raise EngineNotReady("AutoDock Vina is unavailable: install the vina Python package or the Vina binary") from exc
            raise EngineError(f"Vina rescoring failed: {cli_error or exc}") from exc

        try:
            scorer = Vina(sf_name="vina", cpu=1, seed=RANDOM_SEED, verbosity=0)
            scorer.set_receptor(str(receptor_pdbqt))
            scorer.set_ligand_from_file(str(lp))
            scorer.compute_vina_maps(center=list(centroid), box_size=list(size))
            terms = scorer.optimize(max_steps=200)
            energy = float(terms[0])
            scorer.write_pose(str(op), overwrite=True)
            relaxed = op.read_text(encoding="utf-8") if op.exists() else lig_pdbqt
            return energy, relaxed
        except Exception as exc:
            details = f"CLI: {cli_error}; Python API: {exc}" if cli_error else str(exc)
            raise EngineError(f"Vina local rescoring failed: {details}") from exc
