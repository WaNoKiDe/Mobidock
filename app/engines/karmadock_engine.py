# ══ CELL 18 · app/engines/karmadock_engine.py ══
"""
karmadock_engine.py — KarmaDock: graph-neural-network pocket docking (Nat. Comput. Sci. 2023).

Runs the official `utils/virtual_screening.py` from the KarmaDock repo inside its
own Python environment (KarmaDock needs Python 3.9 + PyTorch 1.12 + pyg, which
cannot share an environment with the rest of this app).

Pocket definition: KarmaDock wants a "crystal ligand" MOL2 to locate the pocket.
We write the game's own ligand, moved to the centre of your grid box, as that file.

Energy = Vina score of the relaxed KarmaDock pose (kcal/mol). The KarmaDock
score is kept in the metadata.
"""
import csv
import logging
import math
import os
import subprocess
import tempfile
from pathlib import Path
from typing import List, Optional

from rdkit import Chem

from ..config import ENV_KARMADOCK_DIR, ENV_KARMADOCK_PYTHON, RANDOM_SEED, TOOLS_DIR
from .base import DockingEngine, DockingResult, EngineError
from .poses import first_mol_from_sdf, pose_to_pdb
from .rescoring import vina_relax_and_score, vina_rescoring_available
from .wsl import run_wsl, windows_path_to_wsl, wsl_available

logger = logging.getLogger(__name__)

WEIGHTS = "trained_models/karmadock_screening.pkl"


def karmadock_dir() -> Path:
    return Path(os.environ.get(ENV_KARMADOCK_DIR) or TOOLS_DIR / "KarmaDock")


def karmadock_python() -> Optional[Path]:
    default_env = TOOLS_DIR / "karmadock_env"
    default_py = default_env / "bin" / "python"
    # Do not advertise the half-created environment while its installer is
    # still building PyTorch/ProDy in the startup background thread.
    if not (default_env / ".mobidock-dependencies-ready").exists():
        default_py = ""
    cands = [os.environ.get(ENV_KARMADOCK_PYTHON, ""),
             default_py,
             Path.home() / ".conda/envs/karmadock/bin/python",
             Path("/opt/conda/envs/karmadock/bin/python"),
             Path("/usr/local/envs/karmadock/bin/python")]
    for c in cands:
        if c and Path(c).exists():
            return Path(c)
    return None


def _probe_mol2(ligand_sdf: str, center, path: Path) -> None:
    """
    Write the ligand, translated so its centroid sits at `center`, as a minimal MOL2 file.
    KarmaDock only reads the atom coordinates from it (to locate the pocket).
    """
    mol = Chem.MolFromMolBlock(ligand_sdf, removeHs=False, sanitize=False)
    if mol is None or mol.GetNumAtoms() == 0 or mol.GetNumConformers() == 0:
        raise EngineError("KarmaDock needs a valid 3D ligand structure to locate the binding pocket.")
    conf = mol.GetConformer()
    n = mol.GetNumAtoms()
    cen = [sum(conf.GetAtomPosition(i)[k] for i in range(n)) / n for k in range(3)]
    lines = ["@<TRIPOS>MOLECULE", "pocket_probe", f"{n} 0 0 0 0", "SMALL", "NO_CHARGES", "",
             "@<TRIPOS>ATOM"]
    for i, atom in enumerate(mol.GetAtoms(), start=1):
        p = conf.GetAtomPosition(atom.GetIdx())
        x, y, z = p.x - cen[0] + center[0], p.y - cen[1] + center[1], p.z - cen[2] + center[2]
        sym = atom.GetSymbol()
        lines.append(f"{i:7d} {sym}{i:<5d} {x:10.4f} {y:10.4f} {z:10.4f} {sym:<5s} 1 LIG 0.0000")
    lines += ["@<TRIPOS>BOND", ""]
    Path(path).write_text("\n".join(lines), encoding="utf-8")


import threading
import time


def _clean_env(py, force_cpu: bool = False) -> dict:
    """Isolate KarmaDock's interpreter while preserving a leased server GPU."""
    libdir = str(Path(py).parent.parent / "lib")
    if os.name == "nt":
        env = {"LD_LIBRARY_PATH": windows_path_to_wsl(libdir)}
        if force_cpu:
            env["CUDA_VISIBLE_DEVICES"] = ""
        return env
    env = {k: v for k, v in os.environ.items() if k not in ("PYTHONPATH", "PYTHONHOME", "LD_PRELOAD")}
    env["LD_LIBRARY_PATH"] = libdir + (":" + env["LD_LIBRARY_PATH"] if env.get("LD_LIBRARY_PATH") else "")
    if force_cpu:
        env["CUDA_VISIBLE_DEVICES"] = ""
    return env


def _karmadock_cuda_available(py) -> bool:
    """Check CUDA inside KarmaDock's own interpreter, not just the host process."""
    probe = [str(py), "-c", "import torch; print('MOBIDOCK_CUDA=' + ('1' if torch.cuda.is_available() else '0'))"]
    try:
        if os.name == "nt":
            result = run_wsl(probe, capture_output=True, text=True, timeout=90, env=_clean_env(py))
        else:
            result = subprocess.run(probe, capture_output=True, text=True, timeout=90, env=_clean_env(py))
        return result.returncode == 0 and "MOBIDOCK_CUDA=1" in (result.stdout or "")
    except Exception as exc:
        logger.info("Could not verify KarmaDock CUDA in its interpreter: %s", exc)
        return False


_TORCH = {"ok": None, "msg": "", "ts": 0.0, "running": False}


def _torch_problem() -> List[str]:
    """Checks (in the background, cached 10 min) that `import torch` works inside the KarmaDock env."""
    py = karmadock_python()
    if py is None:
        return []
    st = _TORCH
    if (st["ok"] is None or time.time() - st["ts"] > 600) and not st["running"]:
        st["running"] = True

        def _check():
            try:
                if os.name == "nt":
                    r = run_wsl([str(py), "-c", "import torch"], capture_output=True, text=True,
                                timeout=300, env=_clean_env(py))
                else:
                    r = subprocess.run([str(py), "-c", "import torch"], capture_output=True, text=True,
                                       timeout=300, env=_clean_env(py))
                st["ok"], st["msg"] = (r.returncode == 0), (r.stderr or r.stdout or "")[-300:]
            except Exception as exc:
                st["ok"], st["msg"] = False, str(exc)
            st["ts"], st["running"] = time.time(), False
        threading.Thread(target=_check, daemon=True).start()
    if st["ok"] is False:
        lines = [ln for ln in st["msg"].strip().splitlines() if ln.strip()]
        return ["KarmaDock's PyTorch fails to load (" + (lines[-1][:110] if lines else "unknown error") +
                ") — run Cell 49f (KarmaDock repair), then Cell 49g"]
    return []


class KarmaDockEngine(DockingEngine):
    name = "karmadock"
    display_name = "KarmaDock"
    description = (
        "Deep-learning (graph neural network) pocket-specific docking. "
        "Comparable binding energy = Vina local score of each pose (kcal/mol); KarmaDock's native model score uses different units and is shown separately."
    )
    icon = "⚡"
    color = "#f39c12"
    setup_hint = "On Windows, run the no-WSL Docker shortcut for CPU docking, or use WSL2 for NVIDIA GPU acceleration."

    def problems(self) -> List[str]:
        probs = []
        if os.name == "nt" and not wsl_available():
            probs.append("KarmaDock requires WSL2 with an Ubuntu distribution on Windows")
        d = karmadock_dir()
        if not (d / "utils" / "virtual_screening.py").exists():
            probs.append(f"KarmaDock repo not found at {d} (set KARMADOCK_DIR or run the install cell)")
        elif not (d / WEIGHTS).exists():
            probs.append(f"KarmaDock weights missing: {d / WEIGHTS}")
        if karmadock_python() is None:
            probs.append("KarmaDock Python env not found (set KARMADOCK_PYTHON or run the install cell)")
        if not vina_rescoring_available():
            probs.append("AutoDock Vina is unavailable (install the Vina CLI or Python package to convert KarmaDock poses to kcal/mol)")
        probs.extend(_torch_problem())
        return probs

    def dock(self, receptor_path, ligand_pdbqt, center, size, exhaustiveness=8, **kwargs) -> DockingResult:
        self.require_ready()
        smiles, lig_sdf = kwargs.get("smiles", ""), kwargs.get("ligand_sdf")
        rec_pdb = kwargs.get("receptor_pdb")
        if not smiles or not lig_sdf:
            raise EngineError("KarmaDock needs the ligand SMILES and a 3D ligand.")
        if not rec_pdb or not Path(rec_pdb).exists():
            raise EngineError("Protonated receptor PDB (<name>_H.pdb) not found.")
        if not receptor_path or not Path(receptor_path).exists():
            raise EngineError("Receptor PDBQT is missing (needed for Vina rescoring).")

        # Rebuild a clean isomeric SMILES from the actual 3D input. KarmaDock's
        # own SMILES-to-graph path parses the string repeatedly, and canonical
        # input avoids failures on equivalent but unusual SMILES spellings.
        ligand_mol = Chem.MolFromMolBlock(lig_sdf, removeHs=True, sanitize=True)
        if ligand_mol is None or ligand_mol.GetNumHeavyAtoms() == 0:
            raise EngineError("KarmaDock could not read the 3D ligand as a valid molecule.")
        docking_smiles = Chem.MolToSmiles(ligand_mol, canonical=True, isomericSmiles=True)
        if not docking_smiles:
            raise EngineError("KarmaDock could not create a canonical SMILES for this ligand.")

        repo, py = karmadock_dir(), karmadock_python()
        karmadock_gpu = not kwargs.get("force_cpu") and _karmadock_cuda_available(py)
        n_poses = max(1, min(20, int(kwargs.get("n_poses", 10))))
        with tempfile.TemporaryDirectory() as td:
            td = Path(td)
            protein, smi, mol2, outdir = td / "protein.pdb", td / "ligand.smi", td / "pocket.mol2", td / "out"
            protein.write_text(Path(rec_pdb).read_text(encoding="utf-8"), encoding="utf-8")
            pose_ids = [f"ligand_{index + 1:03d}" for index in range(n_poses)]
            smi.write_text("".join(f"{docking_smiles} {pose_id}\n" for pose_id in pose_ids), encoding="utf-8")
            _probe_mol2(lig_sdf, center, mol2)

            heavy_atoms = max(1, ligand_mol.GetNumHeavyAtoms())
            # Avoid feeding a batch that scales too aggressively for macrocycles
            # and other large ligands. The model is loaded once for all samples.
            atom_budget = max(1, int(os.environ.get("KARMADOCK_BATCH_ATOM_BUDGET", "160")))
            batch_cap = max(1, int(os.environ.get("KARMADOCK_BATCH_SIZE", "4")))
            batch_size = max(1, min(n_poses, batch_cap, atom_budget // heavy_atoms or 1))
            # Tiny ligand graphs are especially sensitive to batch collation in
            # older PyG releases. Run them one-at-a-time; for other inputs retry
            # a failed batch once at size 1 before surfacing the engine error.
            batch_attempts = [1] if heavy_atoms <= 8 else [batch_size]
            if batch_attempts[-1] > 1:
                batch_attempts.append(1)
            res = None
            pose_suffixes = ("_pred_ff_corrected.sdf", "_pred_align_corrected.sdf", "_pred_uncorrected.sdf")
            for batch_index, current_batch in enumerate(batch_attempts):
                if batch_index:
                    for stale in outdir.iterdir():
                        if stale.is_file():
                            stale.unlink()
                cmd = [str(py), "-u", "virtual_screening.py",
                       "--ligand_smi", str(smi), "--protein_file", str(protein),
                       "--crystal_ligand_file", str(mol2), "--out_dir", str(outdir),
                       "--score_threshold=-1000000", "--batch_size", str(current_batch),
                       "--random_seed", str(RANDOM_SEED), "--out_uncoorected", "--out_corrected"]
                if os.name == "nt":
                    res = run_wsl(cmd, cwd=str(repo / "utils"), capture_output=True, text=True,
                                  timeout=1800, env=_clean_env(py, force_cpu=bool(kwargs.get("force_cpu"))))
                else:
                    res = subprocess.run(cmd, cwd=str(repo / "utils"), capture_output=True, text=True,
                                         timeout=1800, env=_clean_env(py, force_cpu=bool(kwargs.get("force_cpu"))))
                has_pose = any((outdir / f"{pose_id}{suffix}").is_file()
                               for pose_id in pose_ids for suffix in pose_suffixes)
                if res.returncode == 0 and has_pose:
                    break
                if batch_index + 1 < len(batch_attempts):
                    logger.warning("KarmaDock returned no usable poses with batch size %s; retrying at batch size 1",
                                   current_batch)
                    continue

            if res is None or res.returncode != 0:
                out = "\n".join(part for part in ((res.stdout if res else "") or "", (res.stderr if res else "") or "") if part.strip())
                if "iJIT_NotifyEvent" in out:
                    raise EngineError("KarmaDock's PyTorch cannot load (MKL version conflict). "
                                      "Run Cell 49f (KarmaDock repair), then Cell 49g.")
                raise EngineError(f"KarmaDock failed: {out[-500:] or 'the docking process returned no details'}")

            scores_by_id = {}
            score_path = outdir / "score.csv"
            if score_path.exists():
                with open(score_path, newline="", encoding="utf-8") as fh:
                    for row in csv.DictReader(fh):
                        scores_by_id[str(row.get("pdb_id") or "")] = row

            records = []
            rescore_errors = []
            for pose_id in pose_ids:
                pose_file = next((outdir / f"{pose_id}{suffix}" for suffix in (
                    "_pred_ff_corrected.sdf", "_pred_align_corrected.sdf", "_pred_uncorrected.sdf"
                ) if (outdir / f"{pose_id}{suffix}").is_file()), None)
                if pose_file is None:
                    continue
                pose = first_mol_from_sdf(pose_file.read_text(encoding="utf-8"))
                if pose is None:
                    continue
                score_row = scores_by_id.get(pose_id, {})
                def score_value(*names):
                    for name in names:
                        try:
                            value = float(score_row.get(name, ""))
                            if math.isfinite(value):
                                return value
                        except (TypeError, ValueError):
                            pass
                    return None
                native_details = {}
                for score_name, raw_value in score_row.items():
                    lowered = str(score_name).lower()
                    if not any(token in lowered for token in ("score", "affinity", "energy")):
                        continue
                    try:
                        parsed_value = float(raw_value)
                    except (TypeError, ValueError):
                        continue
                    if math.isfinite(parsed_value):
                        native_details[score_name] = parsed_value
                try:
                    energy, relaxed = vina_relax_and_score(Path(receptor_path), pose, center, size,
                                                           timeout=int(os.environ.get("KARMADOCK_RESCORING_TIMEOUT", "300")))
                except Exception as exc:
                    logger.warning("Vina rescore failed for %s: %s", pose_id, exc)
                    rescore_errors.append(f"{pose_id}: {exc}")
                    energy, relaxed = None, ""
                records.append({
                    "mol": pose, "file": pose_file.name, "energy": energy, "relaxed": relaxed,
                    "karma_score": score_value("karma_score"),
                    "karma_score_ff": score_value("karma_score_ff"),
                    "karma_score_aligned": score_value("karma_score_aligned"),
                    "native_details": native_details,
                })

            if not records:
                tail = ((res.stdout or "") + (res.stderr or ""))[-500:].strip()
                raise EngineError("KarmaDock produced no pose file. KarmaDock output: " + (tail or "(empty)"))
        scored = [record for record in records if record["energy"] is not None]
        if not scored:
            detail = "; ".join(rescore_errors[-3:])
            suffix = f" Rescore details: {detail}" if detail else ""
            raise EngineError("KarmaDock generated poses, but Vina could not rescore any of them." + suffix)
        records.sort(key=lambda record: (record["energy"] is None, record["energy"] or float("inf")))
        best = next(record for record in records if record["energy"] is not None)
        energy, relaxed = best["energy"], best["relaxed"]

        return DockingResult(
            energy=energy,
            docked_pdbqt=relaxed,
            docked_pdb=pose_to_pdb(relaxed),
            metadata={
                "mode": ("karmadock_gpu" if karmadock_gpu
                         else ("karmadock_cpu_fallback" if kwargs.get("gpu_fallback_reason") else "karmadock_cpu")),
                "gpu_fallback_reason": kwargs.get("gpu_fallback_reason"),
                "karma_score": best["karma_score"],
                "karma_score_ff_corrected": best["karma_score_ff"],
                "karma_score_aligned": best["karma_score_aligned"],
                "native_score_type": ("KarmaDock force-field-corrected score (native units)"
                                      if best["karma_score_ff"] is not None else "KarmaDock score (native units)"),
                "native_details": best["native_details"],
                "pose_source": best["file"],
                "pose_count": len(records), "requested_pose_count": n_poses,
                "energy_source": "Vina local relaxation and score (kcal/mol); same units as AutoDock Vina",
                "native_score_note": "KarmaDock model scores are separate native units and are not kcal/mol.",
                "export_poses": [{
                    "text": Chem.MolToMolBlock(record["mol"]) + "\n$$$$\n",
                    "format": "sdf",
                    "energy": record["karma_score_ff"] if record["karma_score_ff"] is not None else record["karma_score"],
                    "score_type": "KarmaDock score (native units)",
                    "native_score": record["karma_score_ff"] if record["karma_score_ff"] is not None else record["karma_score"],
                    "native_score_type": "KarmaDock force-field-corrected score (native units)" if record["karma_score_ff"] is not None else "KarmaDock score (native units)",
                    "common_be": record["energy"],
                    "common_be_source": "Vina local relaxation and score (kcal/mol)",
                    "common_pose_pdb": pose_to_pdb(record["relaxed"]) if record["relaxed"] else "",
                    "native_details": record["native_details"],
                    "karma_score": record["karma_score"],
                    "karma_score_ff_corrected": record["karma_score_ff"],
                    "karma_score_aligned": record["karma_score_aligned"],
                } for record in records],
            },
            engine_name=self.name,
        )
