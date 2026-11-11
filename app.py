"""Serve the complete MobiDock web app inside a Gradio ZeroGPU Space."""
from __future__ import annotations

import os
import threading
import logging
from pathlib import Path

# ZeroGPU must be imported before MobiDock's optional engine modules.
os.environ["MOBIDOCK_ZERO_GPU"] = "1"
os.environ["MOBIDOCK_PUBLIC_MODE"] = "1"
os.environ["GRADIO_SSR_MODE"] = "false"
# Keep the 1.4 GB GNINA executable in Hugging Face's writable runtime cache;
# the checked-out Space source tree may not be writable by the app process.
_engine_cache = Path(os.environ.get("HF_HOME", "/home/user/.cache/huggingface")) / "mobidock" / "bin"
if not os.environ.get("GNINA_PATH", "").strip():
    os.environ["GNINA_PATH"] = str(_engine_cache / "gnina")

import spaces  # noqa: E402  (required for the ZeroGPU runtime)
import gradio as gr  # noqa: E402

from app.main import app as mobidock_api  # noqa: E402
from app import docking as docking_runtime  # noqa: E402
from app.receptor_prep import RECEPTORS, prepare_all_receptors  # noqa: E402
from app.zero_gpu import register_zero_gpu_runtime, run_engine_in_gpu_lease  # noqa: E402


_BOOTSTRAP_STATE = {"started": False, "finished": False, "errors": []}
_LOGGER = logging.getLogger("mobidock.bootstrap")


def _install_cpu_and_gnina_engines() -> None:
    """Prepare the downloadable engines without delaying the web server startup."""
    try:
        from app.setup_tools import install_gnina, install_karmadock, install_vina

        for name, installer in (("AutoDock Vina", install_vina), ("GNINA", install_gnina), ("KarmaDock", install_karmadock)):
            try:
                attempts = 2 if name == "GNINA" else 1
                for attempt in range(attempts):
                    _LOGGER.info("Setting up %s (attempt %s/%s)", name, attempt + 1, attempts)
                    if installer():
                        break
                    if attempt + 1 < attempts:
                        _LOGGER.warning("%s setup did not complete; retrying", name)
                else:
                    _BOOTSTRAP_STATE["errors"].append(f"{name} installation did not complete")
                    _LOGGER.error("%s installation did not complete", name)
            except Exception as exc:
                _BOOTSTRAP_STATE["errors"].append(f"{name}: {type(exc).__name__}: {exc}")
                _LOGGER.exception("%s installation failed", name)
        # DiffDock is supplied through NVIDIA NIM or an administrator-hosted
        # endpoint; it has no local installer in the laptop/server source.
    finally:
        _BOOTSTRAP_STATE["finished"] = True


def _start_engine_bootstrap() -> None:
    if _BOOTSTRAP_STATE["started"]:
        return
    _BOOTSTRAP_STATE["started"] = True
    if os.environ.get("MOBIDOCK_AUTO_INSTALL_ENGINES", "1").strip().lower() in {"0", "false", "no", "off"}:
        _BOOTSTRAP_STATE["finished"] = True
        return
    threading.Thread(target=_install_cpu_and_gnina_engines, name="mobidock-engine-setup", daemon=True).start()


_start_engine_bootstrap()


def _prepare_public_targets() -> None:
    """Prepare registered targets; this demo starts with none preloaded."""
    try:
        paths = prepare_all_receptors()
        for target_id in RECEPTORS:
            path = paths.get(target_id)
            docking_runtime.RECEPTOR_PATHS[target_id] = path
            if path is None:
                _BOOTSTRAP_STATE["errors"].append(f"{target_id} receptor preparation did not complete")
    except Exception as exc:
        _BOOTSTRAP_STATE["errors"].append(f"Protein setup: {type(exc).__name__}: {exc}")


_prepare_public_targets()

with gr.Blocks(title="MobiDock ZeroGPU runtime") as demo:
    # Both the runtime probe and docking worker are registered Gradio events,
    # because ZeroGPU only leases hardware to event callbacks. FastAPI submits
    # actual engine requests to this hidden docking event via gradio_client.
    _runtime_probe = gr.Button("ZeroGPU runtime check", visible=False)
    _runtime_output = gr.Textbox(visible=False)
    _runtime_probe.click(register_zero_gpu_runtime, inputs=[], outputs=[_runtime_output])

    _gpu_engine_button = gr.Button("Run docking engine on GPU", visible=False)
    _gpu_engine_input = gr.JSON(visible=False)
    _gpu_engine_output = gr.JSON(visible=False)
    _gpu_engine_button.click(
        run_engine_in_gpu_lease,
        inputs=[_gpu_engine_input],
        outputs=[_gpu_engine_output],
        api_name="dock_engine_gpu",
    )

demo.queue(default_concurrency_limit=1, max_size=16)
gradio_host = mobidock_api


if __name__ == "__main__":
    demo.launch(
        _app=gradio_host,
        _frontend=False,
        server_name="0.0.0.0",
        server_port=int(os.environ.get("GRADIO_SERVER_PORT", "7860")),
        ssr_mode=False,
    )
