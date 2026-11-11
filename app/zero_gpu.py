"""ZeroGPU dispatch for docking engines that can use the allocated GPU.

ZeroGPU work must enter through a Gradio event. MobiDock's web UI still uses
its FastAPI job queue, so the API worker submits the engine payload to the
hidden Gradio GPU event exposed by the Space and waits for its result.
"""
from __future__ import annotations

import os
from contextvars import ContextVar
from pathlib import Path
from typing import Any

_GPU_LEASE_ACTIVE: ContextVar[bool] = ContextVar("mobidock_gpu_lease_active", default=False)


def zerogpu_enabled() -> bool:
    return str(
        os.environ.get("MOBIDOCK_ZERO_GPU") or os.environ.get("SPACES_ZERO_GPU") or ""
    ).strip().lower() in {"1", "true", "yes", "on"}


def gpu_lease_active() -> bool:
    """True while the current request is already inside a ZeroGPU event."""
    return _GPU_LEASE_ACTIVE.get()


if zerogpu_enabled():
    # app.py imports spaces before importing any MobiDock engine modules.
    import spaces

    _gpu = spaces.GPU
else:
    # Local laptop/server installations remain independent of Hugging Face.
    def _gpu(function=None, **_options):
        def decorate(fn):
            return fn
        return decorate(function) if function is not None else decorate


def _lease_seconds(request_payload: dict[str, Any]) -> int:
    """Request a short, bounded lease based on the number of requested poses."""
    engine_kwargs = request_payload.get("engine_kwargs") or {}
    try:
        poses = max(1, min(20, int(engine_kwargs.get("n_poses", 10))))
    except (TypeError, ValueError):
        poses = 10
    # Keep individual calls within the short visitor quotas. Larger requests
    # still return as many poses as the engine can produce during this lease.
    return min(120, max(60, 20 + 8 * poses))


def _plain_engine_dock(engine_name: str, engine_kwargs: dict[str, Any]):
    from .engines import get_engine
    return get_engine(engine_name).dock(**engine_kwargs)


def _pathify_engine_kwargs(engine_kwargs: dict[str, Any]) -> dict[str, Any]:
    out = dict(engine_kwargs)
    for key in ("receptor_path", "receptor_pdb"):
        value = out.get(key)
        if value:
            out[key] = Path(value)
    for key in ("center", "size"):
        if key in out:
            out[key] = tuple(out[key])
    return out


def _json_safe(value):
    if isinstance(value, Path):
        return str(value)
    if isinstance(value, dict):
        return {str(key): _json_safe(item) for key, item in value.items()}
    if isinstance(value, (tuple, list)):
        return [_json_safe(item) for item in value]
    return value


@_gpu(duration=_lease_seconds)
def run_engine_in_gpu_lease(request_payload: dict[str, Any]) -> dict[str, Any]:
    """Gradio event handler: run an engine in its allocated ZeroGPU worker."""
    from dataclasses import asdict

    engine_name = str(request_payload.get("engine_name") or "").strip().lower()
    engine_kwargs = request_payload.get("engine_kwargs") or {}
    if engine_name not in {"gnina", "karmadock"} or not isinstance(engine_kwargs, dict):
        raise ValueError("The ZeroGPU docking event received an invalid engine request.")

    token = _GPU_LEASE_ACTIVE.set(True)
    try:
        result = _plain_engine_dock(engine_name, _pathify_engine_kwargs(engine_kwargs))
        return asdict(result)
    finally:
        _GPU_LEASE_ACTIVE.reset(token)


def run_engine_on_gpu(engine_name: str, engine_kwargs: dict[str, Any]):
    """Run directly on local servers, or dispatch through Gradio on ZeroGPU."""
    if not zerogpu_enabled():
        return _plain_engine_dock(engine_name, engine_kwargs)
    if gpu_lease_active():
        return _plain_engine_dock(engine_name, engine_kwargs)

    try:
        from gradio_client import Client
    except ImportError as exc:
        raise RuntimeError("The ZeroGPU docking event client is not installed on this server.") from exc

    port = int(os.environ.get("GRADIO_SERVER_PORT") or os.environ.get("PORT") or "7860")
    client = Client(f"http://127.0.0.1:{port}", verbose=False)
    job = client.submit({"engine_name": engine_name, "engine_kwargs": _json_safe(engine_kwargs)},
                         api_name="/dock_engine_gpu")
    payload = job.result()
    if isinstance(payload, (list, tuple)) and len(payload) == 1:
        payload = payload[0]
    if not isinstance(payload, dict):
        raise RuntimeError("The ZeroGPU docking event returned an invalid engine result.")
    from .engines.base import DockingResult
    return DockingResult(**payload)

@_gpu(duration=1)
def register_zero_gpu_runtime() -> str:
    """A tiny Gradio event that makes the Space's ZeroGPU handler discoverable."""
    return "ZeroGPU docking worker is ready."


def engine_uses_gpu(engine_name: str) -> bool:
    """Return whether this configured backend should request a ZeroGPU lease."""
    if not zerogpu_enabled() or gpu_lease_active():
        return False
    # GNINA's CNN rescoring and KarmaDock's graph model use CUDA when hosted in
    # the Gradio GPU event. Vina stays CPU-based; NIM DiffDock uses NVIDIA's GPU.
    return str(engine_name or "").strip().lower() in {"gnina", "karmadock"}
