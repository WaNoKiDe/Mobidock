"""
routes/synth.py — live synthesizability analysis (used by the designer strip, the Guide and the Dock result card).
"""
import asyncio
import logging
from typing import Optional

from fastapi import APIRouter, HTTPException
from pydantic import BaseModel

from ..synth import SYNTH_MODES, SYNTH_WEIGHT, SA_LIMITS, TIERS, analyze

logger = logging.getLogger("medchem.synth")
router = APIRouter()


class SynthRequest(BaseModel):
    smiles: str
    mode: Optional[str] = "standard"


class SynthExplain3DRequest(BaseModel):
    smiles: str


@router.post("/api/synth/analyze")
async def synth_analyze(req: SynthRequest):
    mode = (req.mode or "standard").strip().lower()
    if mode not in SYNTH_MODES:
        raise HTTPException(status_code=422, detail="mode must be 'standard', 'strict' or 'off'")
    try:
        return await asyncio.to_thread(analyze, (req.smiles or "").strip(), mode)
    except ValueError as exc:
        raise HTTPException(status_code=400, detail=str(exc))
    except Exception as exc:
        logger.exception("Synthesizability analysis failed")
        raise HTTPException(status_code=500, detail=str(exc))


@router.post("/api/synth/explain3d")
async def synth_explain3d(req: SynthExplain3DRequest):
    smiles = (req.smiles or "").strip()
    if not smiles:
        raise HTTPException(status_code=422, detail="SMILES is empty.")
    try:
        from ..sa_explain3d import explain
        return await asyncio.to_thread(explain, smiles)
    except ValueError as exc:
        raise HTTPException(status_code=400, detail=str(exc)) from exc
    except Exception as exc:
        logger.exception("3-D synthesizability explanation failed")
        raise HTTPException(status_code=500, detail="Could not build the synthesizability ligand view.") from exc


@router.get("/api/synth/info")
async def synth_info():
    return {"weight": SYNTH_WEIGHT, "limits": SA_LIMITS,
            "tiers": [{"max_sa": t[0], "key": t[1], "label": t[2], "colour": t[3]} for t in TIERS]}
