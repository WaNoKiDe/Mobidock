# ══ CELL 11 · app/engines/base.py ══
"""
base.py — Abstract base class + shared types for all docking engines.

An engine is "available" only when it can run for real. There are NO mock
results anywhere: if something is missing, `dock()` raises EngineNotReady with
a plain-English list of what to fix, and the UI shows that message.
"""
from abc import ABC, abstractmethod
from dataclasses import dataclass, field
from enum import Enum
from pathlib import Path
from typing import Any, List, Optional, Tuple


class EngineStatus(str, Enum):
    AVAILABLE = "available"        # ready to run for real
    NEEDS_SETUP = "needs_setup"    # something is missing — see `problems`


class EngineError(RuntimeError):
    """The engine ran but failed (bad input, tool crash, API error…)."""


class EngineNotReady(EngineError):
    """The engine cannot run yet — message lists what is missing."""


@dataclass
class DockingResult:
    energy: float                 # kcal/mol (negative = better)
    docked_pdbqt: str             # best pose in the engine's native text format
    metadata: dict[str, Any]
    engine_name: str
    docked_pdb: str = ""          # best pose as PDB text (what the 3D viewer shows)
    is_mock: bool = False         # always False now — kept for API compatibility


class DockingEngine(ABC):
    name: str = "base"
    display_name: str = "Base Engine"
    description: str = ""
    icon: str = "⚙️"
    color: str = "#888"
    setup_hint: str = ""          # one-line "how to enable" shown in the UI tooltip

    # ---- to implement ----------------------------------------------------
    @abstractmethod
    def problems(self) -> List[str]:
        """Return a list of missing requirements ([] means ready)."""

    @abstractmethod
    def dock(
        self,
        receptor_path: Path,
        ligand_pdbqt: str,
        center: Tuple[float, float, float],
        size: Tuple[float, float, float],
        exhaustiveness: int = 8,
        **kwargs,
    ) -> DockingResult:
        """
        kwargs always provided by docking.py:
          smiles        : str
          ligand_sdf    : str   3D MOL block with explicit hydrogens
          receptor_pdb  : Path  protonated receptor PDB (may be None)
        """

    # ---- shared helpers --------------------------------------------------
    def check_availability(self) -> EngineStatus:
        return EngineStatus.NEEDS_SETUP if self.problems() else EngineStatus.AVAILABLE

    @property
    def status(self) -> EngineStatus:
        return self.check_availability()

    def require_ready(self) -> None:
        probs = self.problems()
        if probs:
            raise EngineNotReady(f"{self.display_name} is not set up yet: " + "; ".join(probs))

    def info(self) -> dict:
        probs = self.problems()
        return {
            "name": self.name,
            "display_name": self.display_name,
            "description": self.description,
            "status": (EngineStatus.NEEDS_SETUP if probs else EngineStatus.AVAILABLE).value,
            "problems": probs,
            "setup_hint": self.setup_hint,
            "icon": self.icon,
            "color": self.color,
        }
