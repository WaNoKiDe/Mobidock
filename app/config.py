# ══ CELL 05 · app/config.py ══
"""
config.py — Central paths, environment lookups and constants.

Every other module imports paths from here, so you only ever change a folder
location or an environment-variable name in ONE place.
"""
import os
import shutil
import sys
from pathlib import Path

# ── Folders ───────────────────────────────────────────────────────────────────
BASE_DIR      = Path(__file__).resolve().parent.parent      # project root
APP_DIR       = BASE_DIR / "app"
STATIC_DIR    = BASE_DIR / "static"
DATA_DIR      = BASE_DIR / "data"
RECEPTORS_DIR = BASE_DIR / "receptors"
BIN_DIR       = BASE_DIR / "bin"          # vina / gnina binaries live here
TOOLS_DIR     = BASE_DIR / "tools"        # cloned repos (KarmaDock, DiffDock) live here

for _d in (DATA_DIR, RECEPTORS_DIR, BIN_DIR, TOOLS_DIR):
    _d.mkdir(parents=True, exist_ok=True)

# ── Environment variables you can set to point at your own installs ───────────
ENV_VINA_PATH        = "VINA_PATH"
ENV_GNINA_PATH       = "GNINA_PATH"
ENV_NIM_KEY          = "NVIDIA_NIM_API_KEY"     # DiffDock-L via NVIDIA cloud
ENV_DIFFDOCK_URL     = "DIFFDOCK_URL"           # self-hosted DiffDock NIM endpoint
ENV_DIFFDOCK_DIR     = "DIFFDOCK_DIR"           # local clone of gcorso/DiffDock
ENV_DIFFDOCK_PYTHON  = "DIFFDOCK_PYTHON"
ENV_KARMADOCK_DIR    = "KARMADOCK_DIR"
ENV_KARMADOCK_PYTHON = "KARMADOCK_PYTHON"

# ── Docking defaults ──────────────────────────────────────────────────────────
DEFAULT_EXHAUSTIVENESS = 8
MAX_EXHAUSTIVENESS     = 32
DEFAULT_BOX_SIZE       = (22.0, 22.0, 22.0)
MIN_BOX_EDGE           = 6.0
MAX_BOX_EDGE           = 60.0
RANDOM_SEED            = 42


def has_gpu() -> bool:
    """True if an NVIDIA GPU is visible (nvidia-smi works)."""
    exe = shutil.which("nvidia-smi")
    if not exe:
        return False
    try:
        import subprocess
        r = subprocess.run([exe, "-L"], capture_output=True, text=True, timeout=10)
        return r.returncode == 0 and "GPU" in r.stdout
    except Exception:
        return False


def find_vina() -> "Path | None":
    """Locate the AutoDock Vina CLI binary."""
    cands = [os.environ.get(ENV_VINA_PATH, ""), BIN_DIR / "vina", RECEPTORS_DIR / "vina",
             BIN_DIR / "vina.exe", RECEPTORS_DIR / "vina.exe", shutil.which("vina") or ""]
    for c in cands:
        if c and Path(c).exists():
            return Path(c)
    return None


def python_exe_for(env_var: str) -> str:
    """Python interpreter to use for an external tool (env override or current)."""
    return os.environ.get(env_var) or sys.executable
