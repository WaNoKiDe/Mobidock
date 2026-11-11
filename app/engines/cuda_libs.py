# ══ CELL 14 · app/engines/cuda_libs.py ══
"""
cuda_libs.py — Find CUDA-12 shared libraries so the Gnina binary can start.

The Gnina release binary links against libcudnn.so.9, libcublas.so.12 …
On many GPU installs those libraries ship inside PyTorch's pip packages
(site-packages/nvidia/*/lib) but are NOT on LD_LIBRARY_PATH — that is exactly
the "libcudnn.so.9: cannot open shared object file" error. We add them here.
"""
import glob
import os
import site
import subprocess
import sys
from pathlib import Path
from typing import List


def cuda_lib_dirs() -> List[str]:
    dirs: List[str] = []
    bases = list(dict.fromkeys(site.getsitepackages() + [site.getusersitepackages()] + [p for p in sys.path if p]))
    for base in bases:
        nv = Path(base) / "nvidia"
        if nv.is_dir():
            for sub in sorted(nv.iterdir()):
                lib = sub / "lib"
                if lib.is_dir():
                    dirs.append(str(lib))
        tl = Path(base) / "torch" / "lib"
        if tl.is_dir():
            dirs.append(str(tl))
    from ..config import BIN_DIR
    dirs.append(str(BIN_DIR / "cudalibs"))
    from ..config import DATA_DIR
    wsl_lib_index = DATA_DIR / "gnina_wsl_lib_dirs.txt"
    if wsl_lib_index.is_file():
        try:
            dirs.extend(line.strip() for line in wsl_lib_index.read_text(encoding="utf-8").splitlines() if line.strip())
        except OSError:
            pass
    dirs += glob.glob("/usr/local/cuda*/lib64")
    dirs += [d for d in ("/usr/local/nvidia/lib64", "/usr/lib64-nvidia", "/usr/lib/x86_64-linux-gnu") if os.path.isdir(d)]
    seen, out = set(), []
    for d in dirs:
        if d not in seen:
            seen.add(d)
            out.append(d)
    return out


def env_with_cuda_libs() -> dict:
    env = os.environ.copy()
    dirs = cuda_lib_dirs()
    if os.name == "nt":
        import re
        from .wsl import windows_path_to_wsl
        dirs = [windows_path_to_wsl(d) for d in dirs if not re.match(r"^[A-Za-z]:[/\\]", d)]
    extra = ":".join(dirs)
    suffix = env.get("LD_LIBRARY_PATH", "") if os.name != "nt" else ""
    env["LD_LIBRARY_PATH"] = extra + (":" + suffix if suffix else "")
    return env


def missing_libs(binary: Path) -> List[str]:
    """Names of shared libraries `binary` needs but cannot find (empty list = OK)."""
    try:
        if os.name == "nt":
            from .wsl import run_wsl, windows_path_to_wsl
            r = run_wsl(["ldd", windows_path_to_wsl(binary)], capture_output=True, text=True,
                        env=env_with_cuda_libs(), timeout=30)
        else:
            r = subprocess.run(["ldd", str(binary)], capture_output=True, text=True,
                               env=env_with_cuda_libs(), timeout=30)
    except Exception:
        return []
    return [ln.split("=>")[0].strip() for ln in r.stdout.splitlines() if "not found" in ln]


# library prefix → pip package that provides it
PIP_FOR_LIB = {
    "libcudnn": "nvidia-cudnn-cu12",
    "libcudart": "nvidia-cuda-runtime-cu12",
    "libcublas": "nvidia-cublas-cu12",
    "libcublasLt": "nvidia-cublas-cu12",
    "libcufft": "nvidia-cufft-cu12",
    "libcusparse": "nvidia-cusparse-cu12",
    "libcusolver": "nvidia-cusolver-cu12",
    # CUDA 12.9 removed the legacy NVTX shared object Gnina still links to;
    # pin the compatible CUDA 12.8 wheel for this specific dependency.
    "libnvToolsExt": "nvidia-nvtx-cu12==12.8.55",
    "libnvJitLink": "nvidia-nvjitlink-cu12==12.8.93",
    "libcurand": "nvidia-curand-cu12",
}


def pip_packages_for(missing: List[str]) -> List[str]:
    pk = []
    for lib in missing:
        for prefix, pkg in PIP_FOR_LIB.items():
            if lib.startswith(prefix) and pkg not in pk:
                pk.append(pkg)
    return pk
