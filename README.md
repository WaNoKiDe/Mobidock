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

# MobiDock

MobiDock is a browser-based molecular docking application. This repository contains the Gradio entry point, the FastAPI backend, the web interface, and its Python requirements for a Hugging Face Space.

## Open the hosted app

The current demo is at [wanokide-mobidock-free-demo.hf.space](https://wanokide-mobidock-free-demo.hf.space/).

## Download the Windows website launcher

[Download MobiDock-Web-Setup.exe](downloads/MobiDock-Web-Setup.exe)

The installer adds MobiDock shortcuts to the Windows Desktop and Start menu. Opening a shortcut launches the hosted website in the default browser. It requires an internet connection and does not install or run the docking engines on the computer. This is the website launcher, not the separate local docking installer.

## Create your own Hugging Face Space

### 1. Create the Space

1. Sign in to [Hugging Face](https://huggingface.co/) and open [Create a new Space](https://huggingface.co/new-space).
2. Choose a Space name and owner.
3. Select **Gradio** as the SDK. This repository is not a Docker or Static Space.
4. Choose **Public** if you want anyone to access the app and see or clone its source. Choose a private option if available and you need to restrict access.
5. Choose **ZeroGPU** hardware if it is available to your account. Otherwise, select the hardware available to you. The `suggested_hardware` entry above is only a suggestion; confirm the actual hardware in the Space's **Settings** after creation.

Hugging Face Spaces deploy from Git repositories: a pushed commit triggers a rebuild and restart. See the [Spaces overview](https://huggingface.co/docs/hub/spaces-overview) and [ZeroGPU guide](https://huggingface.co/docs/hub/spaces-zerogpu) for current account eligibility, hardware, and usage limits.

### 2. Copy this repository into the Space

Install [Git for Windows](https://git-scm.com/download/win), then open PowerShell in a folder where you want to keep the two checkouts. Clone both repositories, replacing `HF_USERNAME` and `SPACE_NAME` with the account and Space name you just created:

```powershell
git clone https://github.com/WaNoKiDe/Mobidock.git mobidock-source
git clone https://huggingface.co/spaces/HF_USERNAME/SPACE_NAME.git mobidock-space
robocopy .\mobidock-source .\mobidock-space /E /XD .git downloads __pycache__ /XF *.log
cd .\mobidock-space
git add -A
git commit -m "Deploy MobiDock"
git push
```

`robocopy` may return code `1` when it copied files; that means files were copied. The command excludes Git history, generated Python cache files, logs, and the Windows installer. The Hugging Face clone keeps its own Git remote, so `git push` sends the app files to your Space. When prompted, authenticate with a Hugging Face token that has write access.

The README metadata at the top of this file configures the Space to use Gradio, `app.py`, Gradio 6.31.0, and Python 3.10.13. You can also upload the source files through the Space's **Files** page, but pushing with Git is easier for this many files.

### 3. Configure optional DiffDock access

AutoDock Vina, GNINA, and KarmaDock setup is handled by the app at startup. DiffDock uses an external NVIDIA NIM service or a self-hosted DiffDock endpoint:

- For NVIDIA NIM, add `NVIDIA_NIM_API_KEY` under the Space's **Settings → Variables and secrets → Secrets**.
- For a self-hosted service, set `DIFFDOCK_URL` as a Space variable to that service's endpoint.

These settings are optional for starting the website. Never put API keys in this repository, its README, or a public GitHub commit. Visitors can use their own DiffDock key in the app if the shared Space secret is not configured.

### 4. Wait for the build and check the app

Open the Space's **App** page and wait for the build and startup to finish. The first startup downloads or prepares some docking engines in the background, so their status may take time to become ready. If the build fails, open the Space's **Logs** page and review the build or runtime output.

## Engine and hosting notes

- AutoDock Vina uses CPU docking.
- GNINA and KarmaDock use the Space's ZeroGPU integration when available, with CPU fallback in supported cases.
- DiffDock depends on NVIDIA NIM or an operator-managed DiffDock service; it is not installed as a local engine by this Space.
- Hugging Face hardware availability, queues, and GPU usage limits depend on the account and current platform rules. Free hardware may sleep when unused, and first-time engine setup may need to repeat after a rebuild.
- A public Space exposes its source code and running application to everyone. Keep credentials in Space secrets, not source files.

## Update your Space later

The GitHub repository and the Hugging Face Space are separate Git repositories. Changes pushed to GitHub do not automatically update the Space. To deploy a later GitHub revision, repeat the copy-and-push steps above, or set up a separate GitHub Actions workflow to sync the repositories.
