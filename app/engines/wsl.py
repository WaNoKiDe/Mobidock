"""Small Windows-to-WSL bridge for engines distributed only as Linux programs."""
import os
import re
import shutil
import subprocess
import time
from pathlib import Path
from typing import Optional

_WSL_PROBE = {"checked": 0.0, "available": False}


def wsl_executable() -> Optional[str]:
    return shutil.which("wsl.exe") or shutil.which("wsl")


def windows_path_to_wsl(value) -> str:
    text = str(value).replace("\\", "/")
    match = re.match(r"^([A-Za-z]):/(.*)$", text)
    if match:
        return f"/mnt/{match.group(1).lower()}/{match.group(2)}"
    return text


def wsl_available() -> bool:
    """True when WSL and its default Linux distribution can execute a command."""
    if os.name != "nt":
        return False
    if time.time() - _WSL_PROBE["checked"] < 30:
        return _WSL_PROBE["available"]
    _WSL_PROBE["checked"] = time.time()
    exe = wsl_executable()
    if not exe:
        _WSL_PROBE["available"] = False
        return False
    try:
        result = subprocess.run([exe, "--exec", "true"], capture_output=True, timeout=20)
        _WSL_PROBE["available"] = result.returncode == 0
    except Exception:
        _WSL_PROBE["available"] = False
    return _WSL_PROBE["available"]


def run_wsl(args, *, cwd=None, env=None, user=None, **kwargs):
    """Run a command in WSL, translating Windows workspace paths to /mnt paths."""
    exe = wsl_executable()
    if not exe or not wsl_available():
        raise RuntimeError("WSL2 with a Linux distribution is required for this docking engine.")
    command = [exe]
    if user:
        command.extend(["--user", user])
    if cwd:
        command.extend(["--cd", windows_path_to_wsl(cwd)])
    command.append("--exec")
    if env:
        # Only pass Linux runtime variables. Windows PATH and Python variables do
        # not describe the WSL environment and can break its dynamic loader.
        assignments = []
        for key in ("LD_LIBRARY_PATH", "CUDA_VISIBLE_DEVICES"):
            value = env.get(key)
            if value:
                if key == "LD_LIBRARY_PATH":
                    value = ":".join(windows_path_to_wsl(part) for part in value.split(":"))
                assignments.append(f"{key}={value}")
        if assignments:
            command.extend(["env", *assignments])
    command.extend(windows_path_to_wsl(arg) for arg in args)
    return subprocess.run(command, **kwargs)
