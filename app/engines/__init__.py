# ══ CELL 19 · app/engines/__init__.py ══
"""
engines — Pluggable docking backends.  Registry lives at the bottom of this file.

    vina       AutoDock Vina        (CLI binary or `vina` pip package)
    gnina      Gnina                (CNN-rescored docking, CUDA-12 binary)
    diffdock   DiffDock-L           (NVIDIA NIM cloud API, self-hosted NIM, or local repo)
    karmadock  KarmaDock            (local GNN model in its own Python env)

To add an engine: create <name>_engine.py with a DockingEngine subclass and
add one line to ENGINES below.
"""
from .base import DockingEngine, EngineError, EngineNotReady, EngineStatus
from .vina_engine import VinaEngine
from .gnina_engine import GninaEngine
from .diffdock_engine import DiffDockEngine
from .karmadock_engine import KarmaDockEngine

ENGINES: dict[str, DockingEngine] = {
    "vina":      VinaEngine(),
    "gnina":     GninaEngine(),
    "diffdock":  DiffDockEngine(),
    "karmadock": KarmaDockEngine(),
}


def get_engine(name: str) -> DockingEngine:
    """Get an engine by name. Unknown names raise (no silent fallback to Vina)."""
    key = (name or "vina").lower()
    if key not in ENGINES:
        raise EngineError(f"Unknown docking engine '{name}'. Choose from {list(ENGINES)}")
    return ENGINES[key]


def list_engines() -> list[dict]:
    return [e.info() for e in ENGINES.values()]
