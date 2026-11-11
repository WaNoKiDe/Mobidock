# ══ CELL 30 · app/routes/game.py ══
"""
routes/game.py — Leaderboard, epoch competition, credits ledger, entropy (stress-test) engine, WebSocket.
"""
import asyncio
import logging
import time
from typing import Any, Dict, List, Optional

from fastapi import APIRouter, HTTPException, WebSocket, WebSocketDisconnect

from .. import runtime as rt
from ..runtime import EntropySubmitRequest, SolverSubmission
from ..leaderboard import add_score, get_leaderboard
from ..ledger import (
    append_entry, ResolutionEntry, calculate_reward,
    get_player_credits, get_player_history, get_all_credit_balances,
    deduct_credits, rebuild_credits_from_ledger,
)
from ..entropy_engine import ENTROPY_COST_CREDITS

logger = logging.getLogger("medchem.game")
router = APIRouter()


@router.get("/api/leaderboard")
async def leaderboard(target_id: Optional[str] = None, top_n: int = 10):
    """Return leaderboard scores."""
    return get_leaderboard(target_id=target_id, top_n=top_n)


# ── Routes: Health ────────────────────────────────────────────────────────────



@router.get("/api/v1/epoch/target")
async def get_active_epoch():
    """Fetch the currently active Target Epoch."""
    if rt.epoch_manager is None:
        raise HTTPException(status_code=503, detail="Epoch system not initialised")
    epoch = rt.epoch_manager.get_active_epoch()
    if epoch is None:
        return {"active": False, "epoch": None}
    return {"active": True, "epoch": epoch}


@router.post("/api/v1/epoch/submit")
async def submit_epoch_solution(submission: SolverSubmission):
    """
    Submit a ligand pose for server-side Oracle verification.
    If verified, awards Compute Credits based on solve order.
    """
    if rt.epoch_manager is None or rt.oracle is None:
        raise HTTPException(status_code=503, detail="Systems not initialised")

    # Validate epoch is active
    if not rt.epoch_manager.is_epoch_active(submission.epoch_id):
        raise HTTPException(status_code=410, detail="Epoch expired or invalid")

    current = rt.epoch_manager.current_epoch
    if current is None:
        raise HTTPException(status_code=410, detail="No active epoch")

    solver_id = submission.solver_id.strip() or "Anonymous"

    # Run Oracle verification
    verdict = rt.oracle.verify_pose(
        epoch_id=submission.epoch_id,
        solver_id=solver_id,
        ligand_smiles=current.ligand_smiles,
        ligand_transform=submission.ligand_transform.model_dump(),
        target_id=current.protein_id,
        threshold=current.difficulty_threshold_kcal_mol,
        client_timestamp=submission.client_timestamp or time.time(),
    )

    if not verdict.passed:
        return {
            "status": "rejected",
            "reason": verdict.reason,
            "verified_energy": verdict.verified_energy,
            "collision_count": verdict.collision_count,
            "verification_time_ms": verdict.verification_time_ms,
        }

    # Determine if first solver
    is_first = rt.epoch_manager.mark_solved(submission.epoch_id, solver_id)

    # Calculate time delta from epoch start (or from first solve)
    time_delta = time.time() - current.start_time
    if current.solved_at and not is_first:
        time_delta = time.time() - current.solved_at

    credits = calculate_reward(
        bounty=current.bounty_amount,
        is_first_solve=is_first,
        time_delta_seconds=time_delta,
    )

    # Record in ledger
    entry = ResolutionEntry(
        epoch_id=submission.epoch_id,
        solver_user_id=solver_id,
        docked_pose_hash=verdict.pose_hash,
        verified_binding_energy=verdict.verified_energy,
        credits_awarded=credits,
        timestamp=time.time(),
        is_first_solve=is_first,
        target_id=current.protein_id,
        ligand_smiles=current.ligand_smiles,
    )
    append_entry(entry)

    return {
        "status": "verified",
        "is_first_solve": is_first,
        "credits_awarded": credits,
        "verified_energy": verdict.verified_energy,
        "pose_hash": verdict.pose_hash,
        "collision_count": verdict.collision_count,
        "verification_time_ms": verdict.verification_time_ms,
        "total_credits": get_player_credits(solver_id),
    }


# ── Routes: Player Dashboard ─────────────────────────────────────────────────

@router.get("/api/v1/player/{user_id}/credits")
async def player_credits(user_id: str):
    """Return a player's current Compute Credits balance."""
    return {
        "user_id": user_id,
        "credits": get_player_credits(user_id),
    }


@router.get("/api/v1/player/{user_id}/history")
async def player_history(user_id: str, limit: int = 50):
    """Return a player's resolution history."""
    history = get_player_history(user_id, limit=limit)
    return {
        "user_id": user_id,
        "credits": get_player_credits(user_id),
        "history": history,
    }


@router.get("/api/v1/credits/leaderboard")
async def credits_leaderboard():
    """Return all player credit balances (top earners)."""
    balances = get_all_credit_balances()
    sorted_players = sorted(
        balances.items(), key=lambda x: x[1], reverse=True
    )[:50]
    return {
        "players": [
            {"user_id": uid, "credits": creds}
            for uid, creds in sorted_players
        ]
    }


# ── Routes: Entropy Engine (Phase 4) ─────────────────────────────────────────

@router.post("/api/v1/entropy/submit")
async def submit_entropy_test(req: EntropySubmitRequest):
    """Submit a verified pose for thermodynamic stress-testing."""
    if rt.entropy_engine is None:
        raise HTTPException(status_code=503, detail="Entropy Engine not initialised")

    user_id = req.user_id.strip() or "Anonymous"
    cost = ENTROPY_COST_CREDITS

    # Check credit balance
    current_balance = get_player_credits(user_id)
    if current_balance < cost:
        raise HTTPException(
            status_code=402,
            detail=f"Insufficient credits: {current_balance} < {cost} required",
        )

    # Deduct credits
    if not deduct_credits(user_id, cost):
        raise HTTPException(status_code=402, detail="Credit deduction failed")

    # Record debit in ledger
    debit_entry = ResolutionEntry(
        epoch_id=req.epoch_id,
        solver_user_id=user_id,
        docked_pose_hash=req.pose_hash,
        verified_binding_energy=req.binding_energy,
        credits_awarded=-cost,
        timestamp=time.time(),
        is_first_solve=False,
        target_id=req.target_id,
        ligand_smiles=req.ligand_smiles,
        entry_type="entropy_debit",
    )
    append_entry(debit_entry)

    # Submit job
    job = rt.entropy_engine.submit_job(
        user_id=user_id,
        epoch_id=req.epoch_id,
        target_id=req.target_id,
        pose_hash=req.pose_hash,
        ligand_smiles=req.ligand_smiles,
        binding_energy=req.binding_energy,
        credits_cost=cost,
    )

    # Process asynchronously
    async def _run_and_award():
        completed = await rt.entropy_engine.process_job(job.job_id)
        # If verified, award bonus credits
        if completed.credits_bonus > 0:
            bonus_entry = ResolutionEntry(
                epoch_id=req.epoch_id,
                solver_user_id=user_id,
                docked_pose_hash=req.pose_hash,
                verified_binding_energy=req.binding_energy,
                credits_awarded=completed.credits_bonus,
                timestamp=time.time(),
                is_first_solve=False,
                target_id=req.target_id,
                ligand_smiles=req.ligand_smiles,
                entry_type="entropy_bonus",
            )
            append_entry(bonus_entry)

        # Broadcast result via WebSocket
        dead = []
        for ws in rt.ws_connections:
            try:
                await ws.send_json({
                    "type": "entropy_complete",
                    "data": rt.entropy_engine.get_job(job.job_id),
                })
            except Exception:
                dead.append(ws)
        for ws in dead:
            if ws in rt.ws_connections:
                rt.ws_connections.remove(ws)

    asyncio.create_task(_run_and_award())

    return {
        "status": "queued",
        "job_id": job.job_id,
        "credits_deducted": cost,
        "remaining_credits": get_player_credits(user_id),
    }


@router.get("/api/v1/entropy/job/{job_id}")
async def get_entropy_job(job_id: str):
    """Poll the status of an MD stress-test job."""
    if rt.entropy_engine is None:
        raise HTTPException(status_code=503, detail="Entropy Engine not initialised")
    job = rt.entropy_engine.get_job(job_id)
    if job is None:
        raise HTTPException(status_code=404, detail=f"Job {job_id} not found")
    return job


@router.get("/api/v1/entropy/{user_id}/jobs")
async def get_entropy_jobs(user_id: str, limit: int = 20):
    """List all MD stress-test jobs for a player."""
    if rt.entropy_engine is None:
        raise HTTPException(status_code=503, detail="Entropy Engine not initialised")
    return {
        "user_id": user_id,
        "jobs": rt.entropy_engine.get_player_jobs(user_id, limit=limit),
    }


# ── WebSocket: Real-time Epoch Broadcasts (Phase 2) ──────────────────────────

@router.websocket("/ws/epoch")
async def websocket_epoch(ws: WebSocket):
    """WebSocket endpoint for real-time epoch state updates."""
    await ws.accept()
    rt.ws_connections.append(ws)
    logger.info("WebSocket client connected (%d total)", len(rt.ws_connections))

    # Send current epoch state immediately
    if rt.epoch_manager and rt.epoch_manager.current_epoch:
        try:
            await ws.send_json({
                "type": "current_epoch",
                "data": rt.epoch_manager.get_active_epoch(),
            })
        except Exception:
            pass

    try:
        while True:
            # Keep connection alive — listen for pings
            data = await ws.receive_text()
            if data == "ping":
                await ws.send_json({"type": "pong"})
    except WebSocketDisconnect:
        pass
    finally:
        if ws in rt.ws_connections:
            rt.ws_connections.remove(ws)
        logger.info("WebSocket client disconnected (%d remaining)", len(rt.ws_connections))


# ── Dev runner ────────────────────────────────────────────────────────────────

if __name__ == "__main__":
    import uvicorn
    uvicorn.run(
        "app.main:app",
        host="0.0.0.0",
        port=8000,
        reload=True,
    )
