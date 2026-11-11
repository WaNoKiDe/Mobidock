from ..results_export import sampling_seed, pdbqt_records
# ══ CELL 16 · app/engines/gnina_engine.py ══
"""
gnina_engine.py — Gnina: AutoDock-style docking re-ranked by a 3D CNN.

Runs the real Gnina binary. Reported energy is Gnina's `minimizedAffinity`
(kcal/mol, Vina-scale, so it is comparable with the other engines). The CNN
pose score (0–1) and CNN affinity (pKd) are returned in the metadata.

Works with or without a GPU (adds --no_gpu automatically) but the binary
needs the CUDA-12 shared libraries — see cuda_libs.py.
"""
import logging
import math
import os
import shutil
import subprocess
import tempfile
import time
from pathlib import Path
from typing import List, Optional

from rdkit import Chem

from ..config import BIN_DIR, ENV_GNINA_PATH, RANDOM_SEED, RECEPTORS_DIR, has_gpu
from .base import DockingEngine, DockingResult, EngineError
from .cuda_libs import env_with_cuda_libs, missing_libs, pip_packages_for
from .poses import pose_to_pdb

logger = logging.getLogger(__name__)


def find_gnina() -> Optional[Path]:
    cands = [os.environ.get(ENV_GNINA_PATH, ""), BIN_DIR / "gnina", RECEPTORS_DIR / "gnina", shutil.which("gnina") or ""]
    for c in cands:
        if c and Path(c).is_file():
            return Path(c)
    return None


def _run_gnina(cmd, env, timeout):
    if os.name == "nt" and Path(cmd[0]).suffix.lower() != ".exe":
        from .wsl import run_wsl
        return run_wsl(cmd, capture_output=True, text=True, env=env, timeout=timeout)
    return subprocess.run(cmd, capture_output=True, text=True, env=env, timeout=timeout)


def _is_gpu_failure(output: str) -> bool:
    text = str(output or "").lower()
    markers = (
        "gpu aborted", "gpu abort", "cuda", "cudnn", "cublas", "out of memory",
        "illegal memory access", "device-side assert", "failed to initialize gpu", "aborted",
    )
    return any(marker in text for marker in markers)


def _native_scores(mol):
    scores = {}
    for name in mol.GetPropNames(includePrivate=False, includeComputed=False):
        lowered = name.lower()
        if not any(token in lowered for token in ('score', 'affinity', 'energy')):
            continue
        try:
            value = float(mol.GetProp(name))
        except (TypeError, ValueError):
            continue
        if math.isfinite(value):
            scores[name] = value
    return scores


class GninaEngine(DockingEngine):
    name = "gnina"
    display_name = "Gnina"
    description = (
        "Docking with CNN re-scoring (AutoDock-Vina search + 3D convolutional neural network). "
        "Energy = Gnina minimizedAffinity (kcal/mol); CNN pose score and CNN affinity are "
        "reported alongside."
    )
    icon = "🔮"
    color = "#e74c3c"
    setup_hint = (
        "This public server installs GNINA automatically at startup. The first download is large; "
        "refresh engine status after setup finishes, or restart the Space if it remains unavailable."
        if os.environ.get("MOBIDOCK_PUBLIC_MODE", "").strip().lower() in {"1", "true", "yes", "on"}
        else "Run the 'Install optional engines' cell (downloads Gnina + CUDA-12 libraries)."
    )

    def __init__(self):
        self._cache = (0.0, [])

    def problems(self) -> List[str]:
        if time.time() - self._cache[0] < 30:            # ldd is cheap but UI polls often
            return self._cache[1]
        probs: List[str] = []
        wsl_ok = True
        if os.name == "nt":
            from .wsl import wsl_available
            wsl_ok = wsl_available()
            if not wsl_ok:
                probs.append("Gnina's laptop build runs in WSL2; install WSL2 with an Ubuntu distribution first")
        exe = find_gnina()
        if exe is None:
            if os.environ.get("MOBIDOCK_PUBLIC_MODE", "").strip().lower() in {"1", "true", "yes", "on"}:
                probs.append("GNINA is not installed yet; this server downloads it automatically during startup. Refresh engine status after setup finishes.")
            else:
                probs.append("Gnina binary not found in bin/gnina")
        elif wsl_ok:
            miss = missing_libs(exe)
            if miss:
                pk = " ".join(pip_packages_for(miss))
                probs.append(f"Gnina needs CUDA-12 libraries that are missing: {', '.join(miss)} "
                             f"(fix: pip install --no-deps {pk})")
        self._cache = (time.time(), probs)
        return probs

    def dock(self, receptor_path, ligand_pdbqt, center, size, exhaustiveness=8, **kwargs) -> DockingResult:
        self._cache = (0.0, [])
        self.require_ready()
        exe = find_gnina()
        rec = kwargs.get("receptor_pdb") or receptor_path
        if not rec or not Path(rec).exists():
            raise EngineError("Receptor file is missing — wait for receptor preparation to finish.")
        lig_sdf = kwargs.get("ligand_sdf")
        if not lig_sdf:
            raise EngineError("Gnina needs a 3D ligand (SDF) but none was generated.")

        force_cpu = bool(kwargs.get("force_cpu"))
        gpu = has_gpu() and not force_cpu
        cnn_mode = os.environ.get("GNINA_CNN_SCORING") or ("rescore" if gpu else "none")
        n_poses = max(1, min(20, int(kwargs.get("n_poses", 10))))
        with tempfile.TemporaryDirectory() as td:
            lig, out = Path(td, "ligand.sdf"), Path(td, "docked.sdf")
            lig.write_text(lig_sdf, encoding="utf-8")
            cmd = [
                str(exe), "-r", str(rec), "-l", str(lig), "-o", str(out),
                "--center_x", f"{center[0]:.3f}", "--center_y", f"{center[1]:.3f}", "--center_z", f"{center[2]:.3f}",
                "--size_x", f"{size[0]:.2f}", "--size_y", f"{size[1]:.2f}", "--size_z", f"{size[2]:.2f}",
                "--exhaustiveness", str(exhaustiveness), "--num_modes", str(n_poses),
                "--cnn_scoring", cnn_mode, "--seed", str(sampling_seed.get()),
                "--cpu", str(max(1, min(int(os.environ.get("GNINA_CPU_THREADS", "8")), os.cpu_count() or 2))),
            ]
            if not gpu:
                cmd.append("--no_gpu")
            logger.info("Gnina: %s", " ".join(cmd))
            timeout = int(os.environ.get("GNINA_TIMEOUT", "1800"))
            env = env_with_cuda_libs()
            fallback_reason = str(kwargs.get("gpu_fallback_reason") or "").strip()
            try:
                res = _run_gnina(cmd, env, timeout)
            except subprocess.TimeoutExpired as exc:
                partial = " ".join(str(value or "") for value in (exc.stdout, exc.stderr))
                if not gpu or not _is_gpu_failure(partial):
                    raise EngineError(f"Gnina exceeded its {timeout}s limit: {partial[-300:]}") from exc
                res = None

            if gpu and (res is None or (res.returncode != 0 and _is_gpu_failure((res.stderr or "") + "\n" + (res.stdout or "")))):
                if res is not None:
                    fallback_reason = ((res.stderr or res.stdout or "GNINA GPU process aborted").strip()[-300:])
                logger.warning("GNINA GPU failed; retrying this request in CPU mode: %s", fallback_reason)
                cpu_cmd = [arg for arg in cmd if arg not in ("--no_gpu",)]
                if "--cnn_scoring" in cpu_cmd:
                    score_index = cpu_cmd.index("--cnn_scoring") + 1
                    cpu_cmd[score_index] = "none"
                else:
                    cpu_cmd.extend(["--cnn_scoring", "none"])
                cpu_cmd.extend(["--no_gpu"])
                try:
                    res = _run_gnina(cpu_cmd, env, timeout)
                except subprocess.TimeoutExpired as exc:
                    raise EngineError(f"Gnina GPU aborted and its CPU retry exceeded {timeout}s.") from exc
                gpu = False
                cnn_mode = "none"
            if res.returncode != 0 or not out.exists():
                tail = (res.stderr or res.stdout or "").strip()[-400:]
                raise EngineError(f"Gnina failed (exit {res.returncode}): {tail}")

            sdf_text = out.read_text(encoding="utf-8")
            sup = Chem.SDMolSupplier(str(out), removeHs=False, sanitize=False)
            mols = [m for m in sup if m is not None]
            if not mols:
                raise EngineError("Gnina produced no poses (try a bigger grid box).")
            best = mols[0]                                    # Gnina ranks poses (CNN score under 'rescore')

            def prop(name, default=0.0):
                try:
                    return float(best.GetProp(name))
                except Exception:
                    return default

            affinity = prop("minimizedAffinity")
            first_block = sdf_text.split("$$$$")[0] + "\n$$$$\n"
            pose_pdbs = [Chem.MolToPDBBlock(m) for m in mols]

        return DockingResult(
            energy=affinity,
            docked_pdbqt=first_block,
            docked_pdb=pose_pdbs[0] or pose_to_pdb(first_block),
            metadata={
                "mode": "gnina_gpu" if gpu else ("gnina_cpu_fallback" if fallback_reason else "gnina_cpu"),
                "gpu_fallback_reason": fallback_reason or None,
                "cnn_scoring": cnn_mode,
                "vina_affinity": affinity,
                "cnn_score": prop("CNNscore"),
                "cnn_affinity": prop("CNNaffinity"),
                "n_poses": len(mols), "export_poses": [{
                    "text": Chem.MolToMolBlock(m)+"\n$$$$\n", "format": "sdf",
                    "pdb": pose_pdbs[i],
                    "energy": float(m.GetProp("minimizedAffinity")) if m.HasProp("minimizedAffinity") else None,
                    "score_type": "Gnina minimized affinity (kcal/mol)",
                    "native_score": float(m.GetProp("minimizedAffinity")) if m.HasProp("minimizedAffinity") else None,
                    "native_score_type": "GNINA minimized affinity (kcal/mol)",
                    "native_details": _native_scores(m),
                    "cnn_score": float(m.GetProp("CNNscore")) if m.HasProp("CNNscore") else None,
                    "cnn_affinity": float(m.GetProp("CNNaffinity")) if m.HasProp("CNNaffinity") else None,
                } for i, m in enumerate(mols)],
            },
            engine_name=self.name,
        )
