"""
refstore.py — the game's memory of known docking scores.

Every docking that is run as a yardstick (comparison docking) and every docking a player runs is written here, per
protein, engine, exhaustiveness and grid box. The Ligand Pipeline reads it to answer: "what is the best known docking
score for THIS protein with a molecule like THIS ligand?" A score counts for its protein and engine wherever the molecule was docked (grid box and exhaustiveness are not part of
the lookup); engines are kept apart because their scores are not on the same scale.
"""
import json
import logging
import tempfile
import threading
import time
from pathlib import Path
from typing import Dict, List, Optional

from rdkit import Chem

from .config import DATA_DIR

logger = logging.getLogger("medchem.refstore")
_FILE = DATA_DIR / "reference_scores.json"
_LOCK = threading.RLock()
_DB: Optional[Dict[str, dict]] = None
MAX_RECORDS = 3000


def canon(smiles: str) -> Optional[str]:
    m = Chem.MolFromSmiles((smiles or "").strip())
    return Chem.MolToSmiles(m) if m is not None else None


def box_key(target_id: str) -> str:
    from . import docking
    b = docking.BOX_PARAMS.get(target_id) or {}
    return ",".join(f"{float(v):.1f}" for v in list(b.get("center", ())) + list(b.get("size", ())))


def _load() -> Dict[str, dict]:
    global _DB
    if _DB is None:
        try:
            _DB = json.loads(_FILE.read_text(encoding="utf-8")) if _FILE.exists() else {}
        except Exception as exc:
            logger.warning("could not read %s (%s); starting empty", _FILE, exc)
            _DB = {}
    return _DB


def _save() -> None:
    try:
        _FILE.parent.mkdir(parents=True, exist_ok=True)
        tmp = Path(tempfile.mkstemp(dir=str(_FILE.parent), suffix=".tmp")[1])
        tmp.write_text(json.dumps(_DB), encoding="utf-8")
        tmp.replace(_FILE)
    except Exception as exc:
        logger.warning("could not save %s: %s", _FILE, exc)


def _key(target, cn, engine, exh, box, source, player) -> str:
    return "|".join([target, cn, engine, str(int(exh)), box, source if source == "ref" else f"p:{player}"])


def record(target_id: str, smiles: str, energy: float, engine: str, exhaustiveness: int,
           name: str = "", source: str = "ref", player: str = "", box: str = "") -> Optional[dict]:
    """Stores a docking score. For the same molecule and set-up the BEST (most negative) score is kept, like a record."""
    cn = canon(smiles)
    if cn is None:
        return None
    source = "ref" if source in ("ref", "reference") else ("import" if source == "import" else "player")
    with _LOCK:
        db = _load()
        box = box or box_key(target_id)
        k = _key(target_id, cn, engine, exhaustiveness, box, source, player)
        old = db.get(k)
        rec = {"target_id": target_id, "canon": cn, "smiles": smiles.strip(), "name": name or (old or {}).get("name", ""),
               "energy": float(energy), "engine": engine, "exhaustiveness": int(exhaustiveness), "box": box,
               "source": source, "player": player if source == "player" else "", "runs": 1, "ts": time.time()}
        if old:
            rec["runs"] = old.get("runs", 1) + 1
            rec["energy"] = min(old["energy"], rec["energy"])
        db[k] = rec
        if len(db) > MAX_RECORDS:
            for old_key in sorted(db, key=lambda x: db[x].get("ts", 0))[: len(db) - MAX_RECORDS]:
                db.pop(old_key, None)
        _save()
        return dict(rec)


def query(target_id: str, engine: str, exhaustiveness: int = 0) -> List[dict]:
    """All records for this protein and engine. Grid box and exhaustiveness are NOT filters any more."""
    with _LOCK:
        return [dict(r) for r in _load().values() if r["target_id"] == target_id and r["engine"] == engine]


def get_ref(target_id: str, smiles: str, engine: str, exhaustiveness: int) -> Optional[dict]:
    cn = canon(smiles)
    for r in query(target_id, engine, exhaustiveness):
        if r["source"] == "ref" and r["canon"] == cn:
            return r
    return None


def clear() -> int:
    with _LOCK:
        db = _load()
        n = len(db)
        db.clear()
        _save()
        return n


def best_any(target_id: str, smiles: str, engine: str, sources=None) -> Optional[dict]:
    """The best (most negative) stored score for this molecule on this protein with this engine, wherever it was docked.
    `sources` limits where the score may come from, e.g. ("ref", "import")."""
    cn = canon(smiles)
    rs = [r for r in query(target_id, engine, 0) if r["canon"] == cn and (not sources or r["source"] in sources)]
    return min(rs, key=lambda r: r["energy"]) if rs else None
