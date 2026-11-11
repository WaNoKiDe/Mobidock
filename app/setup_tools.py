# ══ CELL 20 · app/setup_tools.py ══
"""
setup_tools.py — One-call installers for the external docking tools.

Called from the notebook's install cells, e.g.:

    from app.setup_tools import install_vina, install_gnina, install_karmadock
    install_vina()
    install_gnina()
    install_karmadock()

Each function prints what it does and returns True on success.
"""
import os
import shutil
import stat
import subprocess
import sys
import tarfile
import time
import urllib.error
import urllib.request
import zipfile
from pathlib import Path

from .config import BIN_DIR, DATA_DIR, ENV_GNINA_PATH, TOOLS_DIR, has_gpu
from .engines.cuda_libs import env_with_cuda_libs, missing_libs, pip_packages_for
from .engines.wsl import run_wsl, wsl_available

VINA_URL = "https://github.com/ccsb-scripps/AutoDock-Vina/releases/download/v1.2.5/vina_1.2.5_linux_x86_64"
GNINA_BASE = "https://github.com/gnina/gnina/releases/download/v1.3.2/"
GNINA_STD, GNINA_CU128 = "gnina.1.3.2", "gnina.1.3.2.cuda12.8"
KARMADOCK_REPO = "https://github.com/schrojunzhang/KarmaDock.git"
KARMADOCK_ARCHIVE = "https://github.com/schrojunzhang/KarmaDock/archive/refs/heads/main.zip"
MICROMAMBA_URL = "https://micro.mamba.pm/api/micromamba/linux-64/latest"


def _run(cmd, **kw):
    print("  $", " ".join(str(c) for c in cmd))
    if os.name == "nt":
        return run_wsl(cmd, **kw)
    return subprocess.run(cmd, **kw)


def _download(url: str, dest: Path, label: str):
    dest.parent.mkdir(parents=True, exist_ok=True)
    tmp = dest.with_suffix(dest.suffix + ".part")
    last_error = None
    for attempt in range(4):
        offset = tmp.stat().st_size if tmp.exists() else 0
        headers = {"Range": f"bytes={offset}-"} if offset else {}
        request = urllib.request.Request(url, headers=headers)
        try:
            with urllib.request.urlopen(request, timeout=90) as response:
                status = getattr(response, "status", response.getcode())
                content_range = response.headers.get("Content-Range", "")
                # Restart if the host ignored Range; append only when it confirms
                # that the returned bytes start exactly at the saved offset.
                ranged_start = None
                if status == 206 and content_range.startswith("bytes "):
                    try:
                        ranged_start = int(content_range[6:].split("-", 1)[0])
                    except (TypeError, ValueError):
                        ranged_start = None
                append = bool(offset and status == 206 and ranged_start == offset)
                if offset and not append:
                    offset = 0
                mode = "ab" if append else "wb"
                content_length = int(response.headers.get("Content-Length", 0) or 0)
                total = 0
                if content_range and "/" in content_range:
                    try:
                        total = int(content_range.rsplit("/", 1)[1])
                    except ValueError:
                        total = 0
                if not total and content_length:
                    total = offset + content_length
                done, last = offset, -1
                with open(tmp, mode) as file:
                    while True:
                        chunk = response.read(1 << 20)
                        if not chunk:
                            break
                        file.write(chunk)
                        done += len(chunk)
                        if total and int(done * 100 / total) // 10 != last:
                            last = int(done * 100 / total) // 10
                            print(f"    {label}: {done / 1e6:,.0f} / {total / 1e6:,.0f} MB")
                if done <= 0 or (total and done != total):
                    raise IOError(f"incomplete download ({done} of {total or 'unknown'} bytes)")
            os.replace(tmp, dest)
            return
        except urllib.error.HTTPError as exc:
            last_error = exc
            if exc.code == 416 and tmp.exists():
                tmp.unlink(missing_ok=True)
            elif exc.code < 500 and exc.code not in (408, 429):
                raise
        except (urllib.error.URLError, TimeoutError, OSError) as exc:
            last_error = exc
        if attempt < 3:
            print(f"    {label}: download interrupted; resuming ({attempt + 1}/3) …")
            time.sleep(min(2 ** attempt, 8))
    raise IOError(f"Could not finish downloading {label}: {last_error}") from last_error


def _is_gnina_binary(path: Path) -> bool:
    """Reject an HTML error page or a truncated release file before marking setup ready."""
    try:
        if not path.is_file() or path.stat().st_size < 50_000_000:
            return False
        with path.open("rb") as binary:
            return binary.read(4) == b"\x7fELF"
    except OSError:
        return False


def _make_exec(p: Path):
    if os.name == "nt":
        _run(["chmod", "+x", str(p)])
    else:
        p.chmod(p.stat().st_mode | stat.S_IEXEC)


# ── Vina ──────────────────────────────────────────────────────────────────────
def install_vina() -> bool:
    from .config import find_vina
    existing = find_vina()
    if existing:
        print("  AutoDock Vina:", existing.name, "already available")
        return True
    try:
        from vina import Vina  # noqa: F401
        print("  AutoDock Vina: Python API is available for docking and ML-pose rescoring")
        return True
    except Exception:
        pass
    if os.name == "nt":
        print("  AutoDock Vina's Linux installer cannot be used by this Windows Python runtime.")
        return False
    dest = BIN_DIR / "vina"
    if not dest.exists():
        print("Downloading AutoDock Vina 1.2.5 …")
        _download(VINA_URL, dest, "vina")
        _make_exec(dest)
    ok = subprocess.run([str(dest), "--version"], capture_output=True, text=True).returncode == 0
    print("  Vina:", "OK" if ok else "FAILED to run")
    return ok


# ── Gnina ─────────────────────────────────────────────────────────────────────
def _gpu_compute_cap() -> float:
    try:
        out = subprocess.run(["nvidia-smi", "--query-gpu=compute_cap", "--format=csv,noheader"],
                             capture_output=True, text=True, timeout=10).stdout.strip().splitlines()
        return float(out[0]) if out else 0.0
    except Exception:
        return 0.0


def install_gnina() -> bool:
    if os.name == "nt" and not wsl_available():
        print("  Gnina needs WSL2 with an Ubuntu distribution on Windows. Install WSL2, then rerun this cell.")
        return False
    # Public Spaces may mount the source tree read-only. An explicit GNINA_PATH
    # lets deployment keep the large executable in its writable runtime cache.
    dest = Path(os.environ.get(ENV_GNINA_PATH) or (BIN_DIR / "gnina")).expanduser()
    # Hugging Face ZeroGPU does not expose a device to nvidia-smi during
    # startup. Blackwell-class GPUs therefore need the CUDA 12.8 build
    # selected from the deployment mode, not a startup probe.
    zero_gpu_space = str(os.environ.get("MOBIDOCK_ZERO_GPU", "")).strip().lower() in {"1", "true", "yes", "on"}
    asset = GNINA_CU128 if _gpu_compute_cap() >= 10.0 or zero_gpu_space else GNINA_STD
    asset_marker = dest.parent / ".gnina-release-asset"
    installed_asset = asset_marker.read_text(encoding="utf-8").strip() if asset_marker.exists() else ""
    if dest.exists() and not _is_gnina_binary(dest):
        print("  Removing an incomplete or invalid Gnina download …")
        dest.unlink(missing_ok=True)
    if dest.exists() and zero_gpu_space and installed_asset != asset:
        # Replace a previous compatibility build; ZeroGPU startup masks the
        # leased GPU so the normal architecture probe cannot correct it.
        dest.unlink()
    if not dest.exists():
        print(f"Downloading Gnina ({asset}, ~1.4 GB — this takes a few minutes) …")
        _download(GNINA_BASE + asset, dest, "gnina")
        if not _is_gnina_binary(dest):
            dest.unlink(missing_ok=True)
            raise IOError("Downloaded Gnina file is not a complete Linux executable")
        _make_exec(dest)
        asset_marker.write_text(asset + "\n", encoding="utf-8")

    installed: set[str] = set()
    for _ in range(4):
        miss = missing_libs(dest)
        if not miss:
            break
        print("Gnina needs CUDA-12 libraries:", ", ".join(miss))
        pkgs = [p for p in pip_packages_for(miss) if p not in installed]
        if not pkgs:
            break
        installed.update(pkgs)
        print("Installing:", " ".join(pkgs))
        if os.name == "nt":
            libenv = TOOLS_DIR / "gnina_env"
            py = libenv / "bin" / "python"
            if not py.exists():
                created = _run(["python3", "-m", "venv", str(libenv)]).returncode == 0
                if not created:
                    print("Installing Ubuntu's python3-venv support …")
                    run_wsl(["apt-get", "update"], user="root", capture_output=True, text=True, timeout=600)
                    run_wsl(["apt-get", "install", "-y", "python3-venv"], user="root", capture_output=True, text=True, timeout=600)
                    created = _run(["python3", "-m", "venv", str(libenv)]).returncode == 0
                if not created:
                    print("  Could not create the Gnina CUDA-library environment inside Ubuntu.")
                    return False
            install = _run([str(py), "-m", "pip", "install", "-q", "--no-deps", *pkgs])
            if install.returncode == 0:
                libdirs = sorted(str(p) for p in libenv.glob("lib/python*/site-packages/nvidia/*/lib") if p.is_dir())
                if libdirs:
                    (DATA_DIR / "gnina_wsl_lib_dirs.txt").write_text("\n".join(libdirs), encoding="utf-8")
            else:
                break
        else:
            install = _run([sys.executable, "-m", "pip", "install", "-q", "--no-deps", *pkgs])
            if install.returncode != 0:
                break
    miss = missing_libs(dest)
    if miss:
        print("  Gnina: STILL MISSING", miss)
        return False
    r = subprocess.run([str(dest), "--version"], capture_output=True, text=True, env=env_with_cuda_libs())
    print("  Gnina:", (r.stdout or r.stderr).strip().splitlines()[0] if (r.stdout or r.stderr) else "ran")
    print("  GPU detected:", has_gpu(), "(Gnina falls back to CPU automatically)")
    return r.returncode == 0


# ── KarmaDock ─────────────────────────────────────────────────────────────────
def _micromamba() -> Path:
    mm = BIN_DIR / "micromamba"
    if mm.exists():
        return mm
    print("Downloading micromamba …")
    tarball = TOOLS_DIR / "micromamba.tar.bz2"
    _download(MICROMAMBA_URL, tarball, "micromamba")
    with tarfile.open(tarball) as tf:
        member = tf.getmember("bin/micromamba")
        member.name = "micromamba"
        tf.extract(member, BIN_DIR)
    tarball.unlink(missing_ok=True)
    _make_exec(mm)
    return mm


def install_karmadock() -> bool:
    """Clone KarmaDock and build its private Python 3.9 / PyTorch 1.12 environment (~10–15 min)."""
    if os.name == "nt" and not wsl_available():
        print("  KarmaDock needs WSL2 with an Ubuntu distribution on Windows. Install WSL2, then rerun this cell.")
        return False
    repo = TOOLS_DIR / "KarmaDock"
    if not (repo / "utils" / "virtual_screening.py").exists():
        if os.name == "nt":
            print("Downloading KarmaDock source …")
            archive = TOOLS_DIR / "KarmaDock-source.zip"
            _download(KARMADOCK_ARCHIVE, archive, "KarmaDock")
            with zipfile.ZipFile(archive) as source:
                top = next((name.split("/", 1)[0] for name in source.namelist() if "/" in name), "")
                if not top:
                    print("  KarmaDock source archive was empty.")
                    return False
                repo.mkdir(parents=True, exist_ok=True)
                for item in source.infolist():
                    parts = item.filename.split("/")
                    if len(parts) < 2 or parts[0] != top or item.is_dir():
                        continue
                    destination = (repo / Path(*parts[1:])).resolve()
                    if repo.resolve() not in destination.parents:
                        continue
                    destination.parent.mkdir(parents=True, exist_ok=True)
                    with source.open(item) as src, destination.open("wb") as dst:
                        dst.write(src.read())
            archive.unlink(missing_ok=True)
        else:
            print("Cloning KarmaDock …")
            if _run(["git", "clone", "--depth", "1", KARMADOCK_REPO, str(repo)]).returncode != 0:
                return False

    env_dir = TOOLS_DIR / "karmadock_env"
    py = env_dir / "bin" / "python"
    if not py.exists():
        mm = _micromamba()
        # The PyG packages selected by KarmaDock's conda channel are built
        # against CUDA 11.3 even when the dock itself runs on CPU. Always keep
        # the matching runtime libraries so the extensions load; PyTorch still
        # selects CPU execution when no NVIDIA device is available.
        torch_spec = ["pytorch=1.12.1", "cudatoolkit=11.3"]
        mode = "GPU" if has_gpu() else "CPU with CUDA-compatible PyG runtime"
        print(f"Creating KarmaDock environment ({mode}) — be patient …")
        cmd = [str(mm), "create", "-y", "-p", str(env_dir), "-c", "pytorch", "-c", "pyg", "-c", "conda-forge",
               "python=3.9", *torch_spec, "pyg=2.1.0", "mkl<2024.1", "pytorch-scatter", "pytorch-sparse", "pytorch-cluster",
               "pandas", "scipy", "networkx", "tqdm", "joblib"]
        if _run(cmd).returncode != 0:
            print("  micromamba could not solve the environment.")
            return False
    if not py.exists():
        print("  KarmaDock Python environment was not created.")
        return False

    # A Python executable alone does not mean KarmaDock finished installing:
    # the app starts setup in the background, and an early engine-status poll
    # must not race the in-progress Torch/ProDy installation.
    ready_marker = env_dir / ".mobidock-dependencies-ready"
    if not ready_marker.exists():
        pip_pkgs = ["rdkit-pypi==2022.9.1", "prody==2.2.0", "prefetch-generator==1.0.1",
                    "mdanalysis==2.3.0", "rmsd"]
        # ProDy 2.2 needs Cython and NumPy 1.x at build time. PEP 517's
        # isolated build otherwise selects NumPy 2.x and fails to compile.
        if _run([str(py), "-m", "pip", "install", "-q", "numpy<1.24", "Cython<3"]).returncode != 0:
            return False
        if _run([str(py), "-m", "pip", "install", "-q", "--no-build-isolation", *pip_pkgs]).returncode != 0:
            return False
        probe = _run([str(py), "-c", "import torch, torch.torch_version, torch_geometric, prody, MDAnalysis, rdkit"],
                     capture_output=True, text=True)
        if probe.returncode != 0:
            print("  KarmaDock dependency check failed:", (probe.stderr or probe.stdout or "").strip()[-500:])
            return False
        ready_marker.write_text("ready\n", encoding="utf-8")

    os.environ["KARMADOCK_DIR"] = str(repo)
    os.environ["KARMADOCK_PYTHON"] = str(py)
    ok = (repo / "trained_models" / "karmadock_screening.pkl").exists() and py.exists()
    print("  KarmaDock:", "OK" if ok else "incomplete")
    return ok


def status_report() -> None:
    """Print a one-glance status of every engine."""
    from .engines import list_engines
    for e in list_engines():
        mark = "OK" if e["status"] == "available" else "NEEDS SETUP"
        print(f"  [{mark}] {e['display_name']:<14} {e['status']}")
        for p in e["problems"]:
            print(f"        - {p}")
