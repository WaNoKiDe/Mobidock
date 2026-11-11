from ..results_export import dock_with_download
# ══ CELL 27 · app/routes/dock.py ══
"""
routes/dock.py — Docking, engine list/settings, health check.
"""
import asyncio
import threading
import logging
import os
import platform
import sys
import time
from concurrent.futures import ThreadPoolExecutor
from uuid import uuid4

from fastapi import APIRouter, HTTPException
from fastapi.encoders import jsonable_encoder

from .. import docking as dock_engine
from ..config import DATA_DIR, ENV_NIM_KEY, MAX_EXHAUSTIVENESS, find_vina, has_gpu
from ..engines import ENGINES, EngineError, list_engines
from ..leaderboard import add_score
from ..synth import composite_score
from ..runtime import DockRequest, DockResponse, EngineKeyRequest
from ..runtime import DockCompareRequest

logger = logging.getLogger("medchem.dock")
router = APIRouter()
_DOCK_JOB_POOL = ThreadPoolExecutor(
    max_workers=max(1, min(4, int(os.environ.get("MOBIDOCK_DOCK_JOB_WORKERS", "1")))),
    thread_name_prefix="mobidock-dock",
)
_DOCK_JOBS: dict[str, dict] = {}
_DOCK_JOBS_LOCK = threading.Lock()


def _execute_dock_job(job_id: str, req: DockRequest) -> None:
    with _DOCK_JOBS_LOCK:
        job = _DOCK_JOBS.get(job_id)
        if job is None:
            return
        job["status"] = "running"
        job["updated_at"] = time.time()
    try:
        result = asyncio.run(dock_molecule(req))
        payload = jsonable_encoder(result)
        with _DOCK_JOBS_LOCK:
            job = _DOCK_JOBS.get(job_id)
            if job is not None:
                job.update(status="completed", result=payload, updated_at=time.time())
    except HTTPException as exc:
        with _DOCK_JOBS_LOCK:
            job = _DOCK_JOBS.get(job_id)
            if job is not None:
                job.update(status="failed", error=str(exc.detail), updated_at=time.time())
    except Exception as exc:
        logger.exception("Background docking job %s failed", job_id)
        with _DOCK_JOBS_LOCK:
            job = _DOCK_JOBS.get(job_id)
            if job is not None:
                job.update(status="failed", error=f"{type(exc).__name__}: {exc}", updated_at=time.time())


def _prune_dock_jobs() -> None:
    now = time.time()
    for job_id, job in list(_DOCK_JOBS.items()):
        if job.get("status") in {"completed", "failed"} and now - job.get("updated_at", now) > 3600:
            _DOCK_JOBS.pop(job_id, None)


@router.post("/api/dock", response_model=DockResponse)
async def dock_molecule(req: DockRequest):
    """Dock a ligand against a target with the chosen engine (runs in a worker thread)."""
    if req.target_id not in dock_engine.BOX_PARAMS:
        raise HTTPException(status_code=400, detail=f"Unknown target '{req.target_id}'. Select a protein first.")
    smiles = (req.smiles or "").strip()
    if not smiles:
        raise HTTPException(status_code=400, detail="SMILES is empty — load or draw a molecule first.")

    exhaustiveness = max(1, min(req.exhaustiveness or 8, MAX_EXHAUSTIVENESS))
    player = (req.player_name or "").strip() or "Anonymous"
    engine_name = req.engine or "vina"
    lipinski_mode = (req.lipinski_mode or "standard").strip().lower()
    if lipinski_mode not in {"standard", "strict", "off"}:
        raise HTTPException(status_code=422, detail="lipinski_mode must be 'standard', 'strict' or 'off'")
    synth_mode = (req.synth_mode or "standard").strip().lower()
    if synth_mode not in {"standard", "strict", "off"}:
        raise HTTPException(status_code=422, detail="synth_mode must be 'standard', 'strict' or 'off'")

    box = None
    if req.box:
        try:
            center, size = req.box["center"], req.box["size"]
            if len(center) != 3 or len(size) != 3:
                raise ValueError("center and size must each contain three values")
            center, size = dock_engine.clamp_box(center, size)
            box = {"center": center, "size": size}
        except (KeyError, TypeError, ValueError) as exc:
            raise HTTPException(status_code=422, detail=f"Invalid docking box: {exc}") from exc

    try:
        result = await asyncio.to_thread(
            dock_with_download,
            req.target_id, smiles, exhaustiveness, engine_name, lipinski_mode, synth_mode,
            box_override=box, n_poses=max(1, min(20, req.n_poses or 10)),
            diffdock_api_key=req.diffdock_api_key)
    except ValueError as exc:
        raise HTTPException(status_code=422, detail=str(exc))
    except EngineError as exc:
        logger.error("Engine error (%s): %s", engine_name, exc)
        raise HTTPException(status_code=503, detail=str(exc))
    except Exception as exc:
        logger.exception("Docking error")
        raise HTTPException(status_code=500, detail=f"Docking failed: {exc}")

    energy, qed = result["energy"], result["qed"]
    synth = result.get("synth") or {}
    synth_norm = synth.get("norm") if synth.get("available") else None
    composite = composite_score(energy, qed, synth_norm)
    rank = add_score(req.target_id, player, energy, qed, synth_norm, synth.get("sa") if synth.get("available") else None)
    try:
        from .. import refstore
        refstore.record(req.target_id, smiles, energy, (engine_name or "vina").lower(), exhaustiveness, source="player", player=player)
    except Exception as exc:
        import logging as _lg
        _lg.getLogger("medchem.dock").warning("could not store the docking score: %s", exc)
    try:
        from .. import recorddb
        recorddb.save_docking(req.target_id, (engine_name or "vina").lower(), exhaustiveness, player, smiles, energy, qed, composite,
                              result.get("docked_pdb") or "", result.get("props") or {}, synth, result.get("box") or {})
    except Exception as exc:
        import logging as _lg2
        _lg2.getLogger("medchem.dock").warning("could not store the docking in the game database: %s", exc)

    return DockResponse(
        target_id=req.target_id, smiles=smiles, player_name=player,
        energy=energy, qed=qed, score=composite, rank=rank, props=result["props"],
        docked_pdbqt=result["docked_pdbqt"], docked_pdb=result["docked_pdb"],
        ligand_pdbqt=result["ligand_pdbqt"], is_mock=False,
        engine=result["engine"], engine_metadata=result.get("engine_metadata", {}),
        box=result.get("box", {}),
        synth=synth,
    )


@router.post("/api/dock/jobs")
async def submit_docking_job(req: DockRequest):
    """Queue a dock request so long-running ligands don't hold one HTTP connection open."""
    with _DOCK_JOBS_LOCK:
        _prune_dock_jobs()
        if len(_DOCK_JOBS) >= 128:
            raise HTTPException(status_code=429, detail="The docking queue is full. Wait for running jobs to finish and retry.")
        job_id = uuid4().hex
        _DOCK_JOBS[job_id] = {"job_id": job_id, "status": "queued", "created_at": time.time(), "updated_at": time.time()}
    _DOCK_JOB_POOL.submit(_execute_dock_job, job_id, req)
    return {"job_id": job_id, "status": "queued", "poll_url": f"/api/dock/jobs/{job_id}"}


@router.get("/api/dock/jobs/{job_id}")
async def get_docking_job(job_id: str):
    with _DOCK_JOBS_LOCK:
        _prune_dock_jobs()
        job = _DOCK_JOBS.get(job_id)
        if job is None:
            raise HTTPException(status_code=404, detail="Docking job not found; the server may have restarted.")
        return dict(job)


@router.post("/api/dock/rmsd")
async def compare_docked_pose_rmsd(payload: dict):
    """Compare two docked ligand poses in the shared protein coordinate frame."""
    pose_a = str(payload.get("pose_a") or "").strip()
    pose_b = str(payload.get("pose_b") or "").strip()
    if not pose_a or not pose_b:
        raise HTTPException(status_code=422, detail="pose_a and pose_b are required")

    def calculate():
        from ..poses import pose_mol, rmsd_common
        compared = rmsd_common(pose_mol(pose_a), pose_mol(pose_b), min_atoms=1)
        if compared is None:
            return {"available": False, "message": "The poses do not have enough atoms in common for an RMSD comparison."}
        rmsd, atoms = compared
        return {"available": True, "rmsd_angstrom": round(float(rmsd), 3), "atoms_compared": int(atoms)}

    return await asyncio.to_thread(calculate)


@router.post("/api/dock/compare")
async def compare_docking_engines(req: DockCompareRequest):
    """Dock one ligand with each registered engine and export a comparable report."""
    if req.target_id not in dock_engine.BOX_PARAMS:
        raise HTTPException(status_code=400, detail=f"Unknown target '{req.target_id}'. Select a protein first.")
    smiles = (req.smiles or "").strip()
    if not smiles:
        raise HTTPException(status_code=400, detail="SMILES is empty — load or draw a molecule first.")
    lipinski_mode = (req.lipinski_mode or "standard").strip().lower()
    if lipinski_mode not in {"standard", "strict", "off"}:
        raise HTTPException(status_code=422, detail="lipinski_mode must be 'standard', 'strict' or 'off'")
    synth_mode = (req.synth_mode or "standard").strip().lower()
    if synth_mode not in {"standard", "strict", "off"}:
        raise HTTPException(status_code=422, detail="synth_mode must be 'standard', 'strict' or 'off'")
    n_poses = max(1, min(20, int(req.n_poses or 10)))
    exhaustiveness = max(1, min(req.exhaustiveness or 8, MAX_EXHAUSTIVENESS))

    box = None
    if req.box:
        try:
            center, size = dock_engine.clamp_box(req.box["center"], req.box["size"])
            if len(center) != 3 or len(size) != 3:
                raise ValueError("center and size must each contain three values")
            box = {"center": center, "size": size}
        except (KeyError, TypeError, ValueError) as exc:
            raise HTTPException(status_code=422, detail=f"Invalid docking box: {exc}") from exc

    from ..ensemble import run_comparison
    try:
        return await asyncio.to_thread(
            run_comparison,
            req.target_id, smiles, exhaustiveness, lipinski_mode, synth_mode,
            box, n_poses, req.ligand_source, req.ligand_label, req.player_name, req.engine,
            req.diffdock_api_key)
    except ValueError as exc:
        raise HTTPException(status_code=422, detail=str(exc)) from exc
    except Exception as exc:
        logger.exception("Multi-engine docking comparison failed")
        raise HTTPException(status_code=500, detail=f"Docking comparison failed: {exc}") from exc


@router.post("/api/dock/interactions")
async def dock_interactions(payload: dict):
    """Return protein–ligand contacts for drawing in the docked-pose viewer."""
    target_id = str(payload.get("target_id", "")).strip()
    docked_pdb = str(payload.get("docked_pdb", "")).strip()
    if not target_id or not docked_pdb:
        raise HTTPException(status_code=422, detail="target_id and docked_pdb are required")
    from .insight import analyse_interactions
    try:
        return await asyncio.to_thread(analyse_interactions, target_id, docked_pdb)
    except ValueError as exc:
        raise HTTPException(status_code=422, detail=str(exc)) from exc


@router.get("/api/engines")
async def get_engines():
    """Every engine with live status and (if not ready) a list of what is missing."""
    return await asyncio.to_thread(list_engines)


@router.post("/api/engines/diffdock/key")
async def set_diffdock_key(req: EngineKeyRequest):
    """Let the user save an NVIDIA NIM API key for future local server sessions."""
    if os.environ.get("MOBIDOCK_PUBLIC_MODE", "").strip().lower() in {"1", "true", "yes", "on"}:
        raise HTTPException(
            status_code=403,
            detail="This public server accepts its DiffDock key from the server administrator only. Set NVIDIA_NIM_API_KEY in the server environment.",
        )
    key = (req.api_key or "").strip()
    if not key.startswith("nvapi-") or len(key) < 20:
        raise HTTPException(status_code=422, detail="That does not look like an NVIDIA key (should start with 'nvapi-').")
    os.environ[ENV_NIM_KEY] = key
    key_file = DATA_DIR / "nvidia_nim_key.txt"
    key_file.parent.mkdir(parents=True, exist_ok=True)
    key_file.write_text(key, encoding="utf-8")
    return {"status": "ok", "engine": ENGINES["diffdock"].info()}


@router.get("/api/health")
async def health():
    import rdkit
    try:
        import meeko  # noqa: F401
        meeko_ok = True
    except ImportError:
        meeko_ok = False
    try:
        from vina import Vina  # noqa: F401
        vina_py = True
    except Exception:
        vina_py = False
    exe = find_vina()
    engines = await asyncio.to_thread(list_engines)
    return {
        "status": "ok",
        "version": "3.1.0",
        "python": sys.version,
        "platform": platform.platform(),
        "rdkit": rdkit.__version__,
        "vina_python": vina_py,
        "vina_exe": bool(exe),
        "docking_mode": "python_api" if vina_py else ("vina_exe" if exe else "none"),
        "gpu": has_gpu(),
        "meeko": meeko_ok,
        "engines": {e["name"]: e["status"] for e in engines},
        "engine_problems": {e["name"]: e["problems"] for e in engines},
        "receptors": {tid: bool(p and p.exists()) for tid, p in dock_engine.RECEPTOR_PATHS.items()},
    }
