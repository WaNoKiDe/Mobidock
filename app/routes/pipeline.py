
"""
routes/pipeline.py — the Ligand Pipeline.

Every ligand that is loaded goes through, automatically:
    1. route-specific RDKit property warnings for the selected intended use (no drug-likeness score)
    2. synthesizability       (SA score; only if Cell 49k has been run)
    3. closest known molecules (MACCS / Tanimoto against the built-in library)
    4. the BEST KNOWN DOCKING SCORE for this protein and a molecule like this ligand: the number the player has to beat.
       The comparison molecules are docked against the chosen protein with the same engine, grid box and exhaustiveness;
       the scores are remembered (refstore), so the answer also grows with every docking players run.
"""
import asyncio
import inspect
import logging
import threading
from typing import Optional

from fastapi import APIRouter, HTTPException
from pydantic import BaseModel
from rdkit import Chem, DataStructs
from rdkit.Chem import MACCSkeys

from .. import docking as dock_engine
from .. import refstore
from ..engines import EngineError
from ..oral import nearest_known, plan_references
from ..drug_warnings import warnings_all

try:
    from ..synth import analyze as _synth_analyze
except Exception:                       # Cell 49k not installed: the pipeline simply skips this step
    _synth_analyze = None

logger = logging.getLogger("medchem.pipeline")
router = APIRouter()

COMMUNITY_MIN = 0.60                    # other stored molecules count as "similar" from this MACCS Tanimoto upwards
_DOCK_LOCK = threading.Lock()           # one comparison docking at a time, so the player's own docking is not starved


class AnalyzeRequest(BaseModel):
    smiles: str
    target_id: Optional[str] = None


class BestKnownRequest(BaseModel):
    smiles: str
    target_id: str
    engine: Optional[str] = "vina"
    exhaustiveness: int = 8
    player: str = ""


class ReferenceDockRequest(BaseModel):
    target_id: str
    smiles: str
    name: str = "Comparison molecule"
    exhaustiveness: int = 8
    engine: Optional[str] = "vina"
    refresh: bool = False               # True = dock again even if a score is already stored


def _mol(smiles: str):
    mol = Chem.MolFromSmiles((smiles or "").strip())
    if mol is None:
        raise ValueError("That SMILES is not valid")
    return mol


def _analyze(req: AnalyzeRequest) -> dict:
    mol = _mol(req.smiles)
    synth = None
    if _synth_analyze is not None:
        try:
            synth = _synth_analyze(mol, "off")
        except Exception as exc:
            logger.info("synthesizability step skipped: %s", exc)
    prof = warnings_all(mol)
    return {"drug": prof, "synth": synth, "twins": nearest_known(mol, top=5),
            "plan": plan_references(mol, req.target_id)}


_KIND_NOTE = {"same-molecule": "the same molecule", "known-binder": "known binder of this protein",
              "similar-drug": "similar known drug", "community": "docked earlier in this game",
              "extra-reference": "added comparison molecule"}


def _best_known(req: BestKnownRequest) -> dict:
    mol = _mol(req.smiles)
    engine = (req.engine or "vina").strip().lower()
    exh = int(req.exhaustiveness)
    plan = plan_references(mol, req.target_id)
    recs = refstore.query(req.target_id, engine, exh)
    mycanon, me = Chem.MolToSmiles(mol), (req.player or "").strip()
    qfp = MACCSkeys.GenMACCSKeys(mol)

    def usable(cn):
        out = []
        for r in recs:
            if r["canon"] != cn:
                continue
            if r["source"] == "player" and cn == mycanon and r.get("player") == me:
                continue                # you cannot be beaten by your own earlier run of this very molecule
            out.append(r)
        return out

    rows, used = [], set()
    for c in plan["candidates"]:
        used.add(c["canon"])
        rs = usable(c["canon"])
        b = min(rs, key=lambda r: r["energy"]) if rs else None
        rows.append({"name": c["name"], "smiles": c["smiles"], "kind": c["kind"], "kind_note": _KIND_NOTE[c["kind"]],
                     "percent": c["percent"], "note": c["note"], "energy": round(b["energy"], 3) if b else None,
                     "status": "done" if b else "pending", "by": (b or {}).get("player", "")})
    for cn in sorted({r["canon"] for r in recs} - used):
        rs = usable(cn)
        if not rs:
            continue
        b = min(rs, key=lambda r: r["energy"])
        m2 = Chem.MolFromSmiles(b["smiles"])
        sim = DataStructs.TanimotoSimilarity(qfp, MACCSkeys.GenMACCSKeys(m2)) if m2 is not None else 0.0
        if cn != mycanon and sim < COMMUNITY_MIN:
            continue
        kind = "same-molecule" if cn == mycanon else ("community" if b["source"] == "player" else "extra-reference")
        rows.append({"name": b.get("name") or ("docked by " + (b.get("player") or "a player")), "smiles": b["smiles"], "kind": kind,
                     "kind_note": _KIND_NOTE[kind], "percent": round(sim * 100), "note": "", "energy": round(b["energy"], 3),
                     "status": "done", "by": b.get("player", "")})
    done = [r for r in rows if r["energy"] is not None]
    best = min(done, key=lambda r: r["energy"]) if done else None
    return {"target_id": req.target_id, "engine": engine, "exhaustiveness": exh, "rows": rows, "best": best,
            "pending": [{"name": r["name"], "smiles": r["smiles"]} for r in rows if r["status"] == "pending"],
            "loose": plan["loose"], "same": plan["same"], "rationale": plan["rationale"]}


def _reference_dock(req: ReferenceDockRequest) -> dict:
    smiles = (req.smiles or "").strip()
    _mol(smiles)
    engine = (req.engine or "vina").strip().lower()
    exh = int(req.exhaustiveness)
    with _DOCK_LOCK:
        have = refstore.get_ref(req.target_id, smiles, engine, exh)
        if have and not req.refresh:
            return {"name": have["name"] or req.name, "smiles": smiles, "energy": have["energy"], "best_energy": have["energy"],
                    "engine": engine, "exhaustiveness": exh, "cached": True}
        fn = dock_engine.dock_smiles
        params = inspect.signature(fn).parameters
        kwargs = {k: "off" for k in ("lipinski_mode", "synth_mode") if k in params}   # a yardstick is never gated
        result = fn(req.target_id, smiles, exh, engine, **kwargs)
        energy = float(result["energy"])
        rec = refstore.record(req.target_id, smiles, energy, engine, exh, name=req.name, source="ref")
        return {"name": req.name, "smiles": smiles, "energy": energy, "best_energy": rec["energy"] if rec else energy,
                "engine": result.get("engine", engine), "exhaustiveness": exh,
                "seed": (result.get("engine_metadata") or {}).get("seed"), "cached": False}


@router.post("/api/pipeline/analyze")
async def pipeline_analyze(req: AnalyzeRequest):
    try:
        return await asyncio.to_thread(_analyze, req)
    except ValueError as exc:
        raise HTTPException(status_code=422, detail=str(exc))
    except Exception as exc:
        logger.exception("Pipeline analysis failed")
        raise HTTPException(status_code=500, detail=str(exc))


@router.post("/api/pipeline/best_known")
async def pipeline_best_known(req: BestKnownRequest):
    """The best known docking score for this protein and molecules like this ligand (the score to beat)."""
    try:
        return await asyncio.to_thread(_best_known, req)
    except ValueError as exc:
        raise HTTPException(status_code=422, detail=str(exc))
    except Exception as exc:
        logger.exception("best_known failed")
        raise HTTPException(status_code=500, detail=str(exc))


@router.post("/api/pipeline/reference_dock")
async def pipeline_reference_dock(req: ReferenceDockRequest):
    """Docks one comparison molecule against the protein (stored, so repeats are instant)."""
    try:
        return await asyncio.to_thread(_reference_dock, req)
    except ValueError as exc:
        raise HTTPException(status_code=422, detail=str(exc))
    except EngineError as exc:
        raise HTTPException(status_code=503, detail=str(exc))
    except Exception as exc:
        logger.exception("Comparison docking failed")
        raise HTTPException(status_code=500, detail=f"Comparison docking failed: {exc}")


@router.delete("/api/pipeline/reference_cache")
async def pipeline_clear_reference_cache():
    return {"cleared": refstore.clear()}
