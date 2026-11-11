# ══ CELL 17 · app/engines/diffdock_engine.py ══
"""
diffdock_engine.py — DiffDock-L: generative (diffusion) docking.

DiffDock predicts poses with a diffusion model and a confidence score. It does
NOT output kcal/mol, so the chosen pose is locally relaxed + scored with the
Vina scoring function (see rescoring.py). The DiffDock confidence is kept in
the result metadata.

Backends, tried in this order (first one that is configured):
  1. NVIDIA NIM cloud API   → set NVIDIA_NIM_API_KEY  (free key: build.nvidia.com)
  2. Self-hosted DiffDock NIM → set DIFFDOCK_URL (e.g. http://localhost:8000/molecular-docking/diffdock/generate)
  3. Local clone of gcorso/DiffDock → set DIFFDOCK_DIR (+ DIFFDOCK_PYTHON for its env)
"""
import json
import logging
import os
import re
import subprocess
import tempfile
import time
from pathlib import Path
from typing import List, Optional, Tuple

import requests
from rdkit import Chem

from ..config import (ENV_DIFFDOCK_DIR, ENV_DIFFDOCK_PYTHON, ENV_DIFFDOCK_URL, ENV_NIM_KEY,
                      DATA_DIR, has_gpu, python_exe_for)
from .base import DockingEngine, DockingResult, EngineError
from .poses import first_mol_from_sdf, pose_to_pdb
from .rescoring import centroid_in_box, vina_relax_and_score, vina_rescoring_available

logger = logging.getLogger(__name__)

NIM_CLOUD_URL = "https://health.api.nvidia.com/v1/biology/mit/diffdock"
NVCF_STATUS_URL = "https://api.nvcf.nvidia.com/v2/nvcf/pexec/status/{}"


def _backend(api_key: Optional[str] = None) -> Optional[str]:
    if api_key and api_key.strip().startswith("nvapi-") and len(api_key.strip()) >= 20:
        return "nim_cloud"
    if not os.environ.get(ENV_NIM_KEY):
        key_file = DATA_DIR / "nvidia_nim_key.txt"
        if key_file.is_file():
            try:
                key = key_file.read_text(encoding="utf-8").strip()
            except OSError:
                key = ""
            if key.startswith("nvapi-") and len(key) >= 20:
                os.environ[ENV_NIM_KEY] = key
    if os.environ.get(ENV_DIFFDOCK_URL):
        return "nim_selfhosted"
    if os.environ.get(ENV_NIM_KEY):
        return "nim_cloud"
    d = os.environ.get(ENV_DIFFDOCK_DIR)
    if d and (Path(d) / "inference.py").exists():
        return "local"
    return None


def _trim_pdb_around(pdb_text: str, center, radius: float) -> str:
    """Keep only residues with any atom within `radius` Å of center."""
    keep, residues = set(), {}
    for line in pdb_text.splitlines():
        if line.startswith("ATOM"):
            key = (line[21], line[22:27])
            try:
                x, y, z = float(line[30:38]), float(line[38:46]), float(line[46:54])
            except ValueError:
                continue
            residues.setdefault(key, []).append((line, x, y, z))
    out = []
    for key, atoms in residues.items():
        if any((x - center[0]) ** 2 + (y - center[1]) ** 2 + (z - center[2]) ** 2 <= radius ** 2 for _, x, y, z in atoms):
            out.extend(a[0] for a in atoms)
    return ("\n".join(out) + "\nEND\n") if out else (pdb_text or "")


def _pocket_radius(size) -> float:
    """Cover the selected box, plus enough surrounding protein for pocket context."""
    return sum((float(edge) / 2.0 + 8.0) ** 2 for edge in size) ** 0.5


class DiffDockEngine(DockingEngine):
    name = "diffdock"
    display_name = "DiffDock-L"
    description = (
        "Generative diffusion docking (state of the art on PoseBusters). Poses come from DiffDock-L; "
        "energy is the Vina score of the DiffDock pose after local relaxation, DiffDock confidence is "
        "reported alongside."
    )
    icon = "🧬"
    color = "#76b900"
    setup_hint = "Paste a free NVIDIA NIM API key (build.nvidia.com) below, or set DIFFDOCK_URL / DIFFDOCK_DIR."

    def problems(self, api_key: Optional[str] = None) -> List[str]:
        probs = []
        desktop_mode = str(os.environ.get("MOBIDOCK_WINDOWS_DESKTOP", "")).strip().lower() in {"1", "true", "yes", "on"}
        if desktop_mode and os.name == "nt" and not has_gpu():
            probs.append("DiffDock is enabled only on NVIDIA GPU laptops in this Windows build; Vina, GNINA, and KarmaDock remain available on CPU")
        if _backend(api_key) is None:
            probs.append("no backend configured — set NVIDIA_NIM_API_KEY (free key at build.nvidia.com), "
                         "or DIFFDOCK_URL (self-hosted NIM), or DIFFDOCK_DIR (local DiffDock repo)")
        if not vina_rescoring_available():
            probs.append("AutoDock Vina is unavailable (install the Vina CLI or Python package to convert DiffDock poses to kcal/mol)")
        return probs

    # ── main ──────────────────────────────────────────────────────────────────
    def dock(self, receptor_path, ligand_pdbqt, center, size, exhaustiveness=8, **kwargs) -> DockingResult:
        api_key = kwargs.get("diffdock_api_key")
        problems = self.problems(api_key)
        if problems:
            raise EngineError("; ".join(problems))
        smiles = kwargs.get("smiles", "")
        rec_pdb = kwargs.get("receptor_pdb")
        if not smiles:
            raise EngineError("DiffDock needs the ligand SMILES.")
        if not receptor_path or not Path(receptor_path).exists():
            raise EngineError("Receptor PDBQT is missing (needed for Vina rescoring).")
        protein_path = Path(rec_pdb) if rec_pdb else Path(receptor_path).with_name(Path(receptor_path).stem + "_clean.pdb")
        if not protein_path.exists():
            raise EngineError("Receptor PDB file not found for DiffDock.")
        protein_text = protein_path.read_text(encoding="utf-8")

        backend = _backend(api_key)
        t0 = time.time()
        n_poses = max(1, min(20, int(kwargs.get("n_poses", 10))))
        if backend in ("nim_cloud", "nim_selfhosted"):
            poses, confs = self._run_nim(backend, smiles, protein_text, center, size, n_poses, api_key)
        else:
            poses, confs = self._run_local(smiles, protein_path, n_poses, center, size)
        if not poses:
            raise EngineError("DiffDock returned no poses.")

        # Try candidates in the model's confidence order, but retain the Vina
        # local energy for every returned pose. A single malformed top pose
        # must not discard the rest of DiffDock's predictions.
        order = sorted(range(len(poses)), key=lambda i: confs[i] if i < len(confs) else -1e9, reverse=True)
        pose_records = []
        rescore_errors = []
        for i in order:
            in_box = centroid_in_box(poses[i], center, size)
            try:
                energy, relaxed = vina_relax_and_score(
                    Path(receptor_path), poses[i], center, size,
                    timeout=int(os.environ.get("DIFFDOCK_RESCORING_TIMEOUT", "180")),
                )
            except Exception as exc:
                logger.warning("DiffDock pose %d Vina rescore failed: %s", i + 1, exc)
                rescore_errors.append(f"pose {i + 1}: {exc}")
                pose_records.append({"index": i, "energy": None, "relaxed": "", "in_box": in_box,
                                    "rescore_error": str(exc)})
                continue
            pose_records.append({"index": i, "energy": energy, "relaxed": relaxed, "in_box": in_box,
                                "rescore_error": None})

        scored = [record for record in pose_records if record["energy"] is not None]
        if not scored:
            detail = "; ".join(rescore_errors[-3:])
            suffix = f" Rescore details: {detail}" if detail else ""
            raise EngineError("DiffDock returned poses, but Vina could not rescore any of them." + suffix)
        in_box_scored = [record for record in scored if record["in_box"]]
        # DiffDock is a blind-docking model. Prefer poses in the user-selected
        # pocket so a high-scoring pose from another part of the protein cannot
        # become the reported best result.
        best = min(in_box_scored or scored, key=lambda record: record["energy"])
        chosen = best["index"]
        energy, relaxed = best["energy"], best["relaxed"]
        return DockingResult(
            energy=energy,
            docked_pdbqt=relaxed,
            docked_pdb=pose_to_pdb(relaxed),
            metadata={
                "mode": backend,
                "confidence": confs[chosen] if chosen < len(confs) else None,
                "all_confidences": [round(c, 3) for c in confs],
                "pose_count": len(poses), "requested_pose_count": n_poses,
                "export_poses": [{
                    "text": Chem.MolToMolBlock(poses[record["index"]]) + "\n$$$$\n",
                    "pdb": pose_to_pdb(record["relaxed"]) if record["relaxed"] else pose_to_pdb(Chem.MolToMolBlock(poses[record["index"]])),
                    "format": "sdf",
                    "energy": record["energy"] if record["energy"] is not None else (confs[record["index"]] if record["index"] < len(confs) else None),
                    "score_type": "Vina local relaxation and score (kcal/mol)" if record["energy"] is not None else "DiffDock confidence (not kcal/mol)",
                    "common_be": record["energy"],
                    "common_be_source": "Vina local relaxation and score (kcal/mol)",
                    "native_score": confs[record["index"]] if record["index"] < len(confs) else None,
                    "native_score_type": "DiffDock confidence (not kcal/mol)",
                    "confidence": confs[record["index"]] if record["index"] < len(confs) else None,
                    "pose_inside_grid_box": record["in_box"],
                    "rescore_error": record["rescore_error"],
                } for record in pose_records],
                "pose_inside_grid_box": best["in_box"],
                "energy_source": "Vina local relaxation and score of each DiffDock pose",
                "seconds": round(time.time() - t0, 1),
            },
            engine_name=self.name,
        )

    # ── NIM (cloud or self-hosted) ────────────────────────────────────────────
    def _run_nim(self, backend, smiles, protein_text, center, size, n_poses=10,
                 api_key: Optional[str] = None) -> Tuple[List[Chem.Mol], List[float]]:
        url = NIM_CLOUD_URL if backend == "nim_cloud" else (os.environ.get(ENV_DIFFDOCK_URL) or NIM_CLOUD_URL)
        headers = {"Content-Type": "application/json", "Accept": "application/json"}
        key = api_key or os.environ.get(ENV_NIM_KEY)
        if key:
            headers["Authorization"] = f"Bearer {key}"

        def payload(prot: str) -> dict:
            return {
                "ligand": smiles, "ligand_file_type": "txt", "protein": prot,
                "num_poses": n_poses, "time_divisions": 20, "steps": 18,
                "save_trajectory": False, "is_staged": False,
            }

        pocket_text = _trim_pdb_around(protein_text, center, _pocket_radius(size))
        resp = requests.post(url, headers=headers, data=json.dumps(payload(pocket_text)), timeout=300)
        if resp.status_code in (413, 422) and len(protein_text) > 200_000:
            logger.info("DiffDock NIM rejected full protein (%s) — retrying with 30 Å pocket region", resp.status_code)
            trimmed = _trim_pdb_around(protein_text, center, 30.0)
            resp = requests.post(url, headers=headers, data=json.dumps(payload(trimmed)), timeout=300)

        # NVIDIA cloud functions may answer 202 and ask us to poll
        if resp.status_code == 202 and resp.headers.get("NVCF-REQID"):
            req_id = resp.headers["NVCF-REQID"]
            for _ in range(120):
                time.sleep(3)
                resp = requests.get(NVCF_STATUS_URL.format(req_id), headers=headers, timeout=60)
                if resp.status_code != 202:
                    break
        if resp.status_code == 401 or resp.status_code == 403:
            raise EngineError("NVIDIA rejected the API key (401/403). Generate a new key at build.nvidia.com.")
        if resp.status_code != 200:
            raise EngineError(f"DiffDock NIM error {resp.status_code}: {resp.text[:300]}")
        return self._parse_nim(resp.json())

    @staticmethod
    def _parse_nim(data) -> Tuple[List[Chem.Mol], List[float]]:
        raw_poses, confs = [], []
        if isinstance(data, dict):
            raw_poses = data.get("ligand_positions") or data.get("poses") or []
            confs = data.get("position_confidence") or data.get("confidence") or []
        elif isinstance(data, list):                      # older schema: [{ligand_positions, confidence}, …]
            for item in data:
                raw_poses.append(item.get("ligand_positions"))
                confs.append(item.get("confidence", 0.0))
        poses, out_conf = [], []
        for i, rp in enumerate(raw_poses):
            if isinstance(rp, list):
                rp = rp[0] if rp else ""
            mol = first_mol_from_sdf(rp or "")
            if mol is not None:
                poses.append(mol)
                try:
                    out_conf.append(float(confs[i]))
                except Exception:
                    out_conf.append(0.0)
        return poses, out_conf

    # ── Local DiffDock repo ───────────────────────────────────────────────────
    def _run_local(self, smiles, protein_path: Path, n_poses=10, center=None, size=None) -> Tuple[List[Chem.Mol], List[float]]:
        repo = Path(os.environ[ENV_DIFFDOCK_DIR])
        py = python_exe_for(ENV_DIFFDOCK_PYTHON)
        with tempfile.TemporaryDirectory() as td:
            local_protein = protein_path
            if center is not None and size is not None:
                cropped = _trim_pdb_around(
                    protein_path.read_text(encoding="utf-8", errors="ignore"), center, _pocket_radius(size))
                if cropped.strip():
                    local_protein = Path(td, "pocket.pdb")
                    local_protein.write_text(cropped, encoding="utf-8")
            cmd = [py, "-m", "inference", "--config", "default_inference_args.yaml",
                   "--protein_path", str(local_protein), "--ligand_description", smiles,
                   "--complex_name", "game", "--out_dir", td,
                   "--samples_per_complex", str(n_poses), "--batch_size", str(min(5, n_poses))]
            res = subprocess.run(cmd, cwd=str(repo), capture_output=True, text=True, timeout=1800)
            if res.returncode != 0:
                raise EngineError(f"Local DiffDock failed: {(res.stderr or res.stdout)[-400:]}")
            poses, confs = [], []
            for f in sorted(Path(td, "game").glob("rank*_confidence*.sdf")):
                m = re.match(r"rank(\d+)_confidence(-?[\d.]+)\.sdf", f.name)
                mol = first_mol_from_sdf(f.read_text(encoding="utf-8"))
                if m and mol is not None:
                    poses.append(mol)
                    confs.append(float(m.group(2)))
            return poses, confs
