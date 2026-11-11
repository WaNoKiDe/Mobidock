---
title: MobiDock Free Demo
emoji: 🧬
colorFrom: blue
colorTo: purple
sdk: gradio
app_file: app.py
sdk_version: 6.31.0
python_version: "3.10.13"
suggested_hardware: zero-a10g
---

# MobiDock on ZeroGPU

This Space serves the full MobiDock website and its existing APIs: PDB loading, native-ligand selection, molecule editing and viewers, synthesizability and similarity tools, individual engine selection, multi-engine runs, pose tables, RMSD comparisons, reports, ZIP downloads, and Excel exports.

AutoDock Vina runs on CPU. GNINA uses a ZeroGPU GPU lease for docking and CNN rescoring, with CPU fallback if the GPU worker aborts. KarmaDock runs its graph model in a ZeroGPU GPU lease and retries on CPU if the lease fails; it then returns a Vina local-relaxation score in kcal/mol for comparison with AutoDock Vina. KarmaDock's separate native model score uses different units and must not be compared directly with kcal/mol. Longer KarmaDock runs are queued and polled so a large molecule does not have to keep one long browser request open. DiffDock uses NVIDIA NIM or an administrator-hosted DiffDock endpoint. The Space owner can configure `NVIDIA_NIM_API_KEY` as a Space secret or set `DIFFDOCK_URL` as a Space variable. Visitors can enter their own NVIDIA key in the DiffDock engine panel; it is remembered in that browser on their device and sent only with their docking requests.

The first start downloads GNINA and builds KarmaDock's separate CPU environment, so those engine cards may report setup errors until installation finishes. Space storage is temporary and setup may repeat after a container rebuild. DiffDock's external NIM service has its own availability, limits, and key requirements.

ZeroGPU's 24-hour interval is when the daily allowance resets, not a continuous 24-hour GPU session. Hugging Face's current allowance is limited to minutes per day and is shared across visitors. Docking time is also affected by pose count, molecule size, and queue availability. Do not submit private or unpublished molecules to a public Space.
