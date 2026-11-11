# ══ CELL 15 · app/engines/vina_engine.py ══
"""
vina_engine.py — AutoDock Vina (physics-based empirical scoring).

Runs, in order of preference:
  1. the `vina` Python package (pip install vina), if importable
  2. the Vina command-line binary (bin/vina)
If neither works the engine reports exactly what is missing — no fake scores.
"""
import logging
import os
import re
import subprocess
import tempfile
from pathlib import Path
from typing import List, Tuple

from ..config import RANDOM_SEED, find_vina
from .base import DockingEngine, DockingResult, EngineError
from .poses import pose_to_pdb
from ..results_export import sampling_seed

logger = logging.getLogger(__name__)

import threading as _threading

_TL = _threading.local()


def _split_models(text):
    blocks, cur = [], None
    for ln in text.splitlines():
        if ln.startswith("MODEL"):
            cur = [ln]
        elif ln.startswith("ENDMDL") and cur is not None:
            cur.append(ln); blocks.append("\n".join(cur) + "\n"); cur = None
        elif cur is not None:
            cur.append(ln)
    return blocks


def _collect_poses(energies):
    out = []
    for i, b in enumerate(getattr(_TL, "blocks", None) or []):
        try:
            out.append({"energy": energies[i] if i < len(energies) else None, "pdb": pose_to_pdb(b)})
        except Exception:
            pass
    _TL.blocks = []
    return out


def _vina_python():
    try:
        from vina import Vina
        return Vina
    except Exception:
        return None


class VinaEngine(DockingEngine):
    name = "vina"
    display_name = "AutoDock Vina"
    description = (
        "Classical physics-based docking with the AutoDock Vina empirical scoring "
        "function. Fast, reliable and widely validated. Scores are real kcal/mol."
    )
    icon = "⚗️"
    color = "#6c63ff"
    setup_hint = "Run the 'Install' cell — it downloads the Vina binary into bin/."

    def problems(self) -> List[str]:
        if _vina_python() is None and find_vina() is None:
            return ["Vina not found (no `vina` Python package and no bin/vina binary)"]
        return []

    # ── main entry ────────────────────────────────────────────────────────────
    def dock(self, receptor_path, ligand_pdbqt, center, size, exhaustiveness=8, **kwargs) -> DockingResult:
        self.require_ready()
        if not receptor_path or not Path(receptor_path).exists():
            raise EngineError("Receptor PDBQT is missing — wait for receptor preparation to finish.")
        n_poses = max(1, min(20, int(kwargs.get("n_poses", 10) or 10)))

        errors = []
        if _vina_python() is not None:
            try:
                energy, pose, all_e = self._run_python(receptor_path, ligand_pdbqt, center, size, exhaustiveness, n_poses)
                return self._result(energy, pose, "python_api", all_e)
            except Exception as exc:
                logger.warning("Vina Python API failed: %s", exc)
                errors.append(f"python API: {exc}")

        exe = find_vina()
        if exe is not None:
            try:
                energy, pose, all_e = self._run_cli(exe, receptor_path, ligand_pdbqt, center, size, exhaustiveness, n_poses)
                return self._result(energy, pose, "vina_cli", all_e)
            except Exception as exc:
                logger.error("Vina CLI failed: %s", exc)
                errors.append(f"CLI: {exc}")

        raise EngineError("Vina docking failed — " + " | ".join(errors))

    def _result(self, energy, pose, mode, all_energies) -> DockingResult:
        return DockingResult(
            energy=energy, docked_pdbqt=pose, docked_pdb=pose_to_pdb(pose),
            metadata={"mode": mode, "pose_energies": all_energies, "poses": _collect_poses(all_energies)},
            engine_name=self.name,
        )

    # ── Python API ────────────────────────────────────────────────────────────
    @staticmethod
    def _run_python(receptor, ligand_pdbqt, center, size, exh, n_poses) -> Tuple[float, str, list]:
        from vina import Vina
        v = Vina(sf_name="vina", seed=sampling_seed.get(), verbosity=0)
        v.set_receptor(str(receptor))
        with tempfile.TemporaryDirectory() as td:
            lig = Path(td, "ligand.pdbqt")
            out = Path(td, "out.pdbqt")
            lig.write_text(ligand_pdbqt, encoding="utf-8")
            v.set_ligand_from_file(str(lig))
            v.compute_vina_maps(center=list(center), box_size=list(size))
            v.dock(exhaustiveness=exh, n_poses=n_poses)
            # Vina defaults to an energy window of 3 kcal/mol when returning
            # modes, which can silently hide requested poses. Keep all modes
            # found by the search and let the caller report the actual count.
            energy_range = 1000.0
            energies = [round(float(e[0]), 3) for e in v.energies(n_poses=n_poses, energy_range=energy_range)]
            v.write_poses(str(out), n_poses=n_poses, energy_range=energy_range, overwrite=True)
            blocks = _split_models(out.read_text(encoding="utf-8"))
            _TL.blocks = blocks
            return energies[0], (blocks[0] if blocks else out.read_text(encoding="utf-8")), energies

    # ── CLI ───────────────────────────────────────────────────────────────────
    @staticmethod
    def _run_cli(exe, receptor, ligand_pdbqt, center, size, exh, n_poses) -> Tuple[float, str, list]:
        with tempfile.TemporaryDirectory() as td:
            lig, out = Path(td, "ligand.pdbqt"), Path(td, "docked.pdbqt")
            lig.write_text(ligand_pdbqt, encoding="utf-8")
            cmd = [
                str(exe), "--receptor", str(receptor), "--ligand", str(lig), "--out", str(out),
                "--center_x", f"{center[0]:.3f}", "--center_y", f"{center[1]:.3f}", "--center_z", f"{center[2]:.3f}",
                "--size_x", f"{size[0]:.2f}", "--size_y", f"{size[1]:.2f}", "--size_z", f"{size[2]:.2f}",
                "--exhaustiveness", str(exh), "--num_modes", str(n_poses), "--energy_range", "1000", "--seed", str(sampling_seed.get()),
            ]
            res = subprocess.run(cmd, capture_output=True, text=True, timeout=900)
            if res.returncode != 0:
                raise RuntimeError((res.stderr or res.stdout)[-400:].strip())
            energies = []
            for line in res.stdout.splitlines():
                m = re.match(r"\s+(\d+)\s+(-?\d+\.\d+)\s+", line)
                if m:
                    energies.append(float(m.group(2)))
            if not energies:
                raise RuntimeError("could not parse energies from Vina output")
            if not out.exists():
                raise RuntimeError("Vina wrote no output file")
            _txt = out.read_text(encoding="utf-8")
            _TL.blocks = _split_models(_txt)
            return energies[0], _txt, energies
