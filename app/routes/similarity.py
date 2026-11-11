"""
routes/similarity.py — "Drug Twin Finder": how alike are two molecules?

Each molecule gets a fingerprint = a checklist of yes/no chemistry questions (MACCS keys: 166 of them).
    score (Tanimoto) = boxes ticked for BOTH  /  boxes ticked for EITHER      (0 = strangers, 1 = identical)

/api/similarity/compare   rank a list you provide
/api/similarity/discover  AUTO-FIND: your own libraries + PubChem, keep everything above a threshold

Lipinski gates (Rule of Five):
    Gate 1 — BEFORE comparing: molecules that miss too many limits never reach the fingerprints
    Gate 2 — AFTER ranking:    the ranked hits are re-checked (default strict) and failures are dropped
"""
import asyncio
import json
import logging
import time
import urllib.error
import urllib.parse
import urllib.request
from typing import List

from fastapi import APIRouter, HTTPException
from pydantic import BaseModel
from rdkit import Chem, DataStructs
from rdkit.Chem import MACCSkeys

from .. import pubchem
from ..oral import reference_library
from ..lipinski import lipinski_filter
from ..molecule_svg import compare_molecular_features, pair_svgs
from ..synth import synth_summary
from .insight import FEAT_LABEL, FEATS, feature_keys

logger = logging.getLogger("medchem.similarity")
router = APIRouter()

try:
    from rdkit.Chem import rdFingerprintGenerator
    _MORGAN = rdFingerprintGenerator.GetMorganGenerator(radius=2, fpSize=2048)
    def _morgan_fp(m):
        return _MORGAN.GetFingerprint(m)
except Exception:
    from rdkit.Chem import AllChem
    def _morgan_fp(m):
        return AllChem.GetMorganFingerprintAsBitVect(m, 2, nBits=2048)


class Candidate(BaseModel):
    name: str = ""
    smiles: str
    source: str = ""
    cid: int = 0
    chembl: str = ""
    phase: int = -1


class CompareRequest(BaseModel):
    query: str
    candidates: List[Candidate]
    fingerprint: str = "maccs"
    min_similarity: float = 0.0
    top_n: int = 100
    lipinski1_enabled: bool = True
    lipinski1_mode: str = "standard"
    lipinski2_enabled: bool = True
    lipinski2_mode: str = "strict"


class DiscoverRequest(BaseModel):
    query: str
    threshold: float = 0.7
    local: List[Candidate] = []
    use_pubchem: bool = True
    use_reference: bool = True
    use_chembl: bool = True
    max_results: int = 40
    lipinski1_enabled: bool = True
    lipinski1_mode: str = "standard"
    lipinski2_enabled: bool = True
    lipinski2_mode: str = "strict"


def _mol(smiles):
    smiles = (smiles or "").strip()
    return Chem.MolFromSmiles(smiles) if smiles else None


def _level(t: float):
    if t >= 0.90: return "twin", "👯 Twin", 5
    if t >= 0.75: return "sibling", "👫 Sibling", 4
    if t >= 0.55: return "cousin", "🧑‍🤝‍🧑 Cousin", 3
    if t >= 0.35: return "acquaintance", "👋 Acquaintance", 2
    return "stranger", "🧍 Stranger", 1


def _feat(keys):
    return [{"key": k, "emoji": FEAT_LABEL[k][0], "label": FEAT_LABEL[k][1]} for k, *_ in FEATS if k in keys]


def _gate(enabled: bool, mode: str, field: str):
    """Returns (enabled, mode). A gate set to 'off' counts as disabled."""
    mode = (mode or "standard").strip().lower()
    if mode not in ("standard", "strict", "off"):
        raise ValueError(f"{field} must be 'standard', 'strict' or 'off'")
    return (bool(enabled) and mode != "off"), ("standard" if mode == "off" else mode)


def _compare(req: CompareRequest) -> dict:
    g1_on, g1_mode = _gate(req.lipinski1_enabled, req.lipinski1_mode, "lipinski1_mode")
    g2_on, g2_mode = _gate(req.lipinski2_enabled, req.lipinski2_mode, "lipinski2_mode")
    if not 0.0 <= req.min_similarity <= 1.0:
        raise ValueError("min_similarity must be between 0 and 1")
    qm = _mol(req.query)
    if qm is None:
        raise ValueError("Query SMILES is invalid or empty")
    q_lip = lipinski_filter(qm, mode=g1_mode if g1_on else "off")
    if g1_on and not q_lip["passed"]:
        raise ValueError("Lipinski gate 1 rejected your molecule (" + "; ".join(q_lip["violations_detail"]) +
                         "). Untick 'Gate 1' (or choose Standard mode) to compare it anyway.")
    q_canon = Chem.MolToSmiles(qm)
    q_maccs, q_morgan, q_keys = MACCSkeys.GenMACCSKeys(qm), _morgan_fp(qm), feature_keys(qm)
    q_on = set(q_maccs.GetOnBits())
    key = "morgan" if (req.fingerprint or "").strip().lower() == "morgan" else "maccs"

    rows, invalid, rej1 = [], [], []
    for i, c in enumerate(req.candidates[:2500]):
        m = _mol(c.smiles)
        if m is None:
            invalid.append(c.name or c.smiles[:30])
            continue
        name = c.name or f"Compound {i + 1}"
        if g1_on:                                                    # Gate 1: before any fingerprint work
            v1 = lipinski_filter(m, mode=g1_mode)
            if not v1["passed"]:
                rej1.append({"index": i, "name": name, "smiles": c.smiles.strip(), "reasons": v1["violations_detail"]})
                continue
        c_maccs, c_morgan, c_keys = MACCSkeys.GenMACCSKeys(m), _morgan_fp(m), feature_keys(m)
        c_on = set(c_maccs.GetOnBits())
        t_maccs = DataStructs.TanimotoSimilarity(q_maccs, c_maccs)
        t_morgan = DataStructs.TanimotoSimilarity(q_morgan, c_morgan)
        same = Chem.MolToSmiles(m) == q_canon
        score = t_morgan if key == "morgan" else t_maccs
        lvl, lvl_label, stars = _level(score)
        rows.append({
            "index": i, "name": name, "smiles": c.smiles.strip(),
            "source": c.source, "cid": c.cid, "chembl": c.chembl, "phase": c.phase, "same": same,
            "maccs": round(t_maccs, 8), "morgan": round(t_morgan, 8),
            "percent": round(score * 100, 8), "fingerprint_identical": score == 1.0,
            "level": lvl, "level_label": lvl_label, "stars": stars,
            "shared_features": _feat(q_keys & c_keys), "only_query": _feat(q_keys - c_keys),
            "only_candidate": _feat(c_keys - q_keys),
            "bits_query": len(q_on), "bits_candidate": len(c_on), "bits_shared": len(q_on & c_on),
            "ro5": lipinski_filter(m, mode="standard"),
            "synth": synth_summary(m),
        })

    rows = [r for r in rows if r[key] >= req.min_similarity]
    rows.sort(key=lambda r: r[key], reverse=True)

    rej2 = []
    if g2_on:                                                        # Gate 2: after similarity ranking
        kept = []
        for r in rows:
            v2 = lipinski_filter(_mol(r["smiles"]), mode=g2_mode)
            if v2["passed"]:
                kept.append(r)
            else:
                rej2.append({"index": r["index"], "name": r["name"], "smiles": r["smiles"], "percent": r["percent"],
                             "reasons": v2["violations_detail"]})
        rows = kept

    rows = rows[:max(1, req.top_n)]
    for rank, r in enumerate(rows, 1):
        candidate = _mol(r["smiles"])
        try:
            r["query_svg"], r["candidate_svg"] = pair_svgs(qm, candidate)
            r["group_features"] = compare_molecular_features(qm, candidate)
        except Exception as exc:
            logger.debug("Could not draw similarity structures for %s: %s", r.get("name"), exc)
            r["query_svg"], r["candidate_svg"] = "", ""
            r["group_features"] = {"shared": [], "only_query": [], "only_candidate": []}
        r["rank"] = rank
    return {"query": {"smiles": req.query.strip(), "canonical": q_canon, "maccs_bits": len(q_on),
                      "morgan_bits": q_morgan.GetNumOnBits(), "features": _feat(q_keys), "ro5": lipinski_filter(qm, "standard")},
            "ranked_by": key, "results": rows, "invalid": invalid,
            "lipinski": {"gate1": {"enabled": g1_on, "mode": g1_mode, "rejected": rej1},
                         "gate2": {"enabled": g2_on, "mode": g2_mode, "rejected": rej2}}}


# ─────────────────────────── online sources: PubChem + ChEMBL (robust, with real error messages) ───────────────────────────
class NetError(Exception):
    """A failed web request, with a message a player can act on."""
    code = 0


_RETRY_CODES = {429, 500, 502, 503, 504}
_CHEMBL = "https://www.ebi.ac.uk/chembl/api/data"


def _http_json(url: str, data: bytes = None, timeout: int = 25, tries: int = 3, label: str = "request"):
    headers = {**pubchem.UA, "Accept": "application/json"}
    if data is not None:
        headers["Content-Type"] = "application/x-www-form-urlencoded"
    last = None
    for attempt in range(tries):
        try:
            req = urllib.request.Request(url, data=data, headers=headers)
            with urllib.request.urlopen(req, timeout=timeout) as r:
                body = r.read().decode("utf-8", "ignore").strip()
                return json.loads(body) if body else {}
        except urllib.error.HTTPError as exc:
            try:
                detail = exc.read().decode("utf-8", "ignore").replace("\n", " ")[:140]
            except Exception:
                detail = ""
            last = NetError(f"{label}: HTTP {exc.code} {exc.reason}" + (f" — {detail}" if detail else ""))
            last.code = exc.code
            if exc.code in _RETRY_CODES and attempt < tries - 1:        # server busy: wait and try again
                time.sleep(1.5 * (attempt + 1))
                continue
            raise last
        except (urllib.error.URLError, TimeoutError, OSError, ValueError) as exc:
            last = NetError(f"{label}: {type(exc).__name__} — {getattr(exc, 'reason', exc)}")
            if attempt < tries - 1:
                time.sleep(1.0 * (attempt + 1))
                continue
            raise last
    raise last


def _pc_cids(d: dict) -> list:
    return list((d.get("IdentifierList") or {}).get("CID", []))


def _pc_strategies(smiles: str, threshold: int, max_records: int):
    """Different ways to ask PubChem the same question. SMILES with / or \\ cannot go in a URL path, so they use POST."""
    q = urllib.parse.quote(smiles, safe="")
    opts = f"Threshold={int(threshold)}&MaxRecords={int(max_records)}"
    body = urllib.parse.urlencode({"smiles": smiles}).encode()
    path_ok = "/" not in smiles and "\\" not in smiles
    out = []
    if path_ok:
        out.append(("fastsimilarity_2d (GET)", f"{pubchem.BASE}/fastsimilarity_2d/smiles/{q}/cids/JSON?{opts}", None))
    out.append(("fastsimilarity_2d (POST)", f"{pubchem.BASE}/fastsimilarity_2d/smiles/cids/JSON?{opts}", body))
    if path_ok:
        out.append(("similarity (GET)", f"{pubchem.BASE}/similarity/smiles/{q}/cids/JSON?{opts}", None))
    out.append(("similarity (POST)", f"{pubchem.BASE}/similarity/smiles/cids/JSON?{opts}", body))
    return out


def _pc_follow(d: dict, wait: int = 40) -> list:
    """Some PubChem searches run as a background job: poll the list key until the CIDs are ready."""
    key = (d.get("Waiting") or {}).get("ListKey")
    t0 = time.time()
    while key and time.time() - t0 < wait:
        time.sleep(1.5)
        d = _http_json(f"{pubchem.BASE}/listkey/{key}/cids/JSON", timeout=30, tries=2, label="PubChem (waiting for job)")
        if "IdentifierList" in d:
            break
        key = (d.get("Waiting") or {}).get("ListKey") or key
    return _pc_cids(d)


def pubchem_similar_cids(smiles: str, pc_threshold: int, max_records: int = 150):
    """PubChem compounds similar to `smiles` (PubChem's own 2-D fingerprint). Returns (cids, method, failures)."""
    failures = []
    for label, url, body in _pc_strategies(smiles, pc_threshold, max_records):
        try:
            return _pc_follow(_http_json(url, body, timeout=30, tries=2, label=f"PubChem {label}")), label, failures
        except NetError as exc:
            if exc.code == 404:                       # "no matches" is an answer, not a failure
                return [], label, failures
            failures.append(str(exc))
    try:                                              # last resort: look the molecule up by CID, then search by CID
        d = _http_json(f"{pubchem.BASE}/smiles/cids/JSON", urllib.parse.urlencode({"smiles": smiles}).encode(),
                       timeout=25, tries=2, label="PubChem CID lookup")
        cid = (_pc_cids(d) or [0])[0]
        if cid:
            url = f"{pubchem.BASE}/fastsimilarity_2d/cid/{cid}/cids/JSON?Threshold={int(pc_threshold)}&MaxRecords={int(max_records)}"
            return _pc_follow(_http_json(url, timeout=30, tries=2, label="PubChem by CID")), "by CID", failures
    except NetError as exc:
        if exc.code != 404:
            failures.append(str(exc))
    raise NetError(f"all {len(_pc_strategies(smiles, 0, 0)) + 1} PubChem methods failed. Last error: " +
                   (failures[-1] if failures else "unknown"))


def pubchem_props(cids: list) -> list:
    out = []
    for i in range(0, len(cids), 100):
        chunk = ",".join(str(c) for c in cids[i:i + 100])
        for props in ("Title,SMILES,ConnectivitySMILES", "Title,IsomericSMILES,CanonicalSMILES"):
            try:
                d = _http_json(f"{pubchem.BASE}/cid/{chunk}/property/{props}/JSON", timeout=30, tries=2, label="PubChem properties")
                out += d.get("PropertyTable", {}).get("Properties", [])
                break
            except NetError as exc:
                if exc.code != 400:
                    raise
    return out


def chembl_similar(smiles: str, percent: int, limit: int = 40) -> list:
    """Compounds ChEMBL finds similar (ECFP4 Tanimoto >= percent). Includes the drug-development phase (4 = approved)."""
    url = f"{_CHEMBL}/similarity/{urllib.parse.quote(smiles, safe='')}/{max(40, min(100, int(percent)))}.json?limit={int(limit)}"
    d = _http_json(url, timeout=45, tries=2, label="ChEMBL")
    out = []
    for m in d.get("molecules", []) or []:
        smi = (m.get("molecule_structures") or {}).get("canonical_smiles")
        if not smi:
            continue
        try:
            phase = int(float(m.get("max_phase"))) if m.get("max_phase") not in (None, "") else -1
        except (TypeError, ValueError):
            phase = -1
        out.append({"name": m.get("pref_name") or m.get("molecule_chembl_id") or "ChEMBL compound", "smiles": smi,
                    "chembl": m.get("molecule_chembl_id") or "", "phase": phase})
    return out


def _discover(req: DiscoverRequest) -> dict:
    qm = _mol(req.query)
    if qm is None:
        raise ValueError("Your molecule's SMILES is invalid or empty")
    th = max(0.3, min(0.99, req.threshold))
    notes, pool = [], {}
    counts = {"library": 0, "known": 0, "pubchem": 0, "chembl": 0}

    def add(name, smiles, source, bucket, cid=0, chembl="", phase=-1):
        m = _mol(smiles)
        if m is None:
            return
        canon = Chem.MolToSmiles(m)
        if canon not in pool:
            pool[canon] = Candidate(name=name, smiles=smiles, source=source, cid=cid, chembl=chembl, phase=phase)
            counts[bucket] += 1

    for c in req.local:
        add(c.name or c.smiles[:20], c.smiles, c.source or "library", "library")
    if req.use_reference:                                  # always available, works offline
        for name, smi in reference_library():
            add(name, smi, "known-drug", "known")
    if req.use_pubchem:
        try:
            cids, how, _ = pubchem_similar_cids(req.query.strip(), max(50, round(th * 100) - 15), 150)
            for p in pubchem_props(cids):
                smi = pubchem.pick_smiles(p)
                if smi:
                    add(p.get("Title") or f"CID {p.get('CID')}", smi, "pubchem", "pubchem", cid=int(p.get("CID") or 0))
            notes.append(f"🌐 PubChem: OK via {how} — {len(cids)} compounds proposed.")
        except Exception as exc:
            logger.warning("PubChem similarity failed: %s", exc)
            notes.append(f"⚠️ PubChem: {exc}")
    if req.use_chembl:
        try:
            found = chembl_similar(req.query.strip(), max(40, round(th * 100) - 10), 40)
            for f in found:
                add(f["name"], f["smiles"], "chembl", "chembl", chembl=f["chembl"], phase=f["phase"])
            notes.append(f"🌐 ChEMBL: OK — {len(found)} compounds proposed.")
        except Exception as exc:
            logger.warning("ChEMBL similarity failed: %s", exc)
            notes.append(f"⚠️ ChEMBL: {exc}")
    if (req.use_pubchem or req.use_chembl) and not counts["pubchem"] and not counts["chembl"]:
        notes.append("Showing known drugs from the built-in library only. Press 🔌 Test connection to see exactly what is blocked.")

    # ChEMBL's similarity endpoint uses ECFP4/Morgan. Re-rank and apply the
    # selected threshold with the same fingerprint so valid ChEMBL hits are
    # not discarded by a second, unrelated MACCS threshold.
    chembl_only = req.use_chembl and not req.use_pubchem and not req.use_reference and not req.local
    result = _compare(CompareRequest(
        query=req.query, candidates=list(pool.values()),
        fingerprint="morgan" if chembl_only else "maccs", min_similarity=th, top_n=req.max_results + 5,
        lipinski1_enabled=req.lipinski1_enabled, lipinski1_mode=req.lipinski1_mode,
        lipinski2_enabled=req.lipinski2_enabled, lipinski2_mode=req.lipinski2_mode))
    result["results"] = [r for r in result["results"] if not r["same"]][:req.max_results]
    for rank, r in enumerate(result["results"], 1):
        r["rank"] = rank
    result.update({"threshold": th, "notes": notes, "n_checked": len(pool), "n_library": counts["library"] + counts["known"],
                   "n_pubchem": counts["pubchem"], "n_chembl": counts["chembl"]})
    return result


def _diagnose() -> dict:
    """Runs each way of reaching the online databases once and reports exactly what happened."""
    smi = "CC(=O)Oc1ccccc1C(=O)O"          # aspirin
    rows = []

    def attempt(label, fn):
        t0 = time.time()
        try:
            rows.append({"name": label, "ok": True, "detail": fn(), "ms": round((time.time() - t0) * 1000)})
        except Exception as exc:
            rows.append({"name": label, "ok": False, "detail": str(exc)[:240], "ms": round((time.time() - t0) * 1000)})

    def pc_basic():
        d = _http_json(f"{pubchem.BASE}/cid/2244/property/MolecularFormula/JSON", timeout=15, tries=1, label="PubChem")
        return "reachable — " + d["PropertyTable"]["Properties"][0]["MolecularFormula"]

    def make(url, body):
        def run():
            d = _http_json(url, body, timeout=25, tries=1, label="PubChem")
            if "Waiting" in d:
                return "accepted (runs as a background job)"
            return f"{len(_pc_cids(d))} compounds"
        return run

    attempt("PubChem · basic lookup", pc_basic)
    for label, url, body in _pc_strategies(smi, 90, 5):
        attempt(f"PubChem · {label}", make(url, body))
    attempt("ChEMBL · similarity", lambda: f"{len(chembl_similar(smi, 90, 5))} compounds")
    ok = [r for r in rows if r["ok"]]
    return {"checks": rows, "pubchem_ok": any(r["ok"] and r["name"].startswith("PubChem · ") and "basic" not in r["name"] for r in rows),
            "chembl_ok": any(r["ok"] and r["name"].startswith("ChEMBL") for r in rows), "any_ok": bool(ok)}


@router.post("/api/similarity/compare")
async def compare_similarity(req: CompareRequest):
    if not req.candidates:
        raise HTTPException(status_code=400, detail="Provide at least one candidate molecule")
    try:
        return await asyncio.to_thread(_compare, req)
    except ValueError as exc:
        raise HTTPException(status_code=400, detail=str(exc))
    except Exception as exc:
        logger.exception("Similarity failed")
        raise HTTPException(status_code=500, detail=str(exc))


@router.post("/api/similarity/discover")
async def discover_similar(req: DiscoverRequest):
    """AUTO-FIND every molecule at or above the look-alike threshold (own libraries + PubChem)."""
    try:
        return await asyncio.to_thread(_discover, req)
    except ValueError as exc:
        raise HTTPException(status_code=400, detail=str(exc))
    except Exception as exc:
        logger.exception("Discovery failed")
        raise HTTPException(status_code=500, detail=str(exc))


@router.get("/api/similarity/diagnose")
async def diagnose_online():
    """Tests every route to PubChem and ChEMBL and says what worked (used by the 🔌 Test connection button)."""
    return await asyncio.to_thread(_diagnose)
