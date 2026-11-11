"""
recorddb.py — the game's database (SQLite file data/medchem.db).

  molecules   every molecule players have docked: properties, drug score, synthesizability, closest known drug, parent
  dockings    every docking a player ran: score, pose, box, quality numbers
  poses       the pose of each reference (blind) docking, for the RMSD comparison
  records     the chain of "score to beat" records: each time a player beats the reference, their molecule is stored here
  claims      the verdict for each docking (so asking twice gives the same answer)
  natives     the ligand pose of the experimental PDB entry (used to judge poses when two energies tie)
"""
import json
import sqlite3
import threading
from contextlib import closing
from datetime import datetime
from typing import Any, Dict, List, Optional

from rdkit import Chem

from .config import DATA_DIR

DB_FILE = DATA_DIR / "medchem.db"
_LOCK = threading.RLock()
_SCHEMA = """
CREATE TABLE IF NOT EXISTS molecules(
  canon TEXT PRIMARY KEY, smiles TEXT, first_player TEXT, first_target TEXT, first_seen TEXT, times_docked INTEGER DEFAULT 0,
  props TEXT, synth TEXT, qed REAL, drug_score REAL, drug_grade TEXT, nearest_name TEXT, nearest_pct INTEGER,
  known_drug INTEGER DEFAULT 0, parent_canon TEXT, parent_smiles TEXT);
CREATE TABLE IF NOT EXISTS dockings(
  id INTEGER PRIMARY KEY AUTOINCREMENT, target_id TEXT, engine TEXT, exhaustiveness INTEGER, player TEXT, canon TEXT, smiles TEXT,
  energy REAL, qed REAL, score REAL, pose_pdb TEXT, box TEXT, created TEXT);
CREATE INDEX IF NOT EXISTS ix_dockings ON dockings(target_id, engine, canon);
CREATE TABLE IF NOT EXISTS poses(
  target_id TEXT, canon TEXT, engine TEXT, energy REAL, pose_pdb TEXT, source TEXT, created TEXT, PRIMARY KEY(target_id, canon, engine));
CREATE TABLE IF NOT EXISTS records(
  id INTEGER PRIMARY KEY AUTOINCREMENT, target_id TEXT, engine TEXT, parent_canon TEXT, holder_canon TEXT, holder_smiles TEXT,
  holder_name TEXT, player TEXT, energy REAL, prev_energy REAL, delta REAL, rmsd_to_prev REAL, decided_by TEXT, pose_pdb TEXT,
  docking_id INTEGER, created TEXT);
CREATE INDEX IF NOT EXISTS ix_records ON records(target_id, engine, parent_canon, id);
CREATE TABLE IF NOT EXISTS claims(docking_id INTEGER PRIMARY KEY, verdict TEXT, created TEXT);
CREATE TABLE IF NOT EXISTS natives(target_id TEXT PRIMARY KEY, comp TEXT, pdb TEXT, created TEXT);
"""


def _now() -> str:
    return datetime.now().strftime("%Y-%m-%d %H:%M:%S")


def _conn() -> sqlite3.Connection:
    DB_FILE.parent.mkdir(parents=True, exist_ok=True)
    c = sqlite3.connect(str(DB_FILE), timeout=30)
    c.row_factory = sqlite3.Row
    c.executescript(_SCHEMA)
    return c


def _row(r) -> Optional[Dict[str, Any]]:
    return dict(r) if r is not None else None


def canon(smiles: str) -> str:
    m = Chem.MolFromSmiles(smiles or "")
    return Chem.MolToSmiles(m) if m is not None else ""


# ── dockings + molecules ─────────────────────────────────────────────────────
def save_docking(target_id, engine, exhaustiveness, player, smiles, energy, qed, score, pose_pdb, props, synth, box) -> int:
    """Store one player docking and register the molecule (with all its properties) if it is new."""
    cn = canon(smiles)
    prof = nearest = None
    try:
        from .oral import drug_profile, nearest_known
        m = Chem.MolFromSmiles(smiles)
        prof = drug_profile(m)
        nearest = (nearest_known(m, top=1) or [None])[0]
    except Exception:
        pass
    with _LOCK, closing(_conn()) as c:
        cur = c.execute("INSERT INTO dockings(target_id,engine,exhaustiveness,player,canon,smiles,energy,qed,score,pose_pdb,box,created)"
                        " VALUES(?,?,?,?,?,?,?,?,?,?,?,?)",
                        (target_id, engine, int(exhaustiveness), player, cn, smiles, float(energy), qed, score, pose_pdb or "",
                         json.dumps(box or {}), _now()))
        did = cur.lastrowid
        if c.execute("SELECT 1 FROM molecules WHERE canon=?", (cn,)).fetchone():
            c.execute("UPDATE molecules SET times_docked=times_docked+1 WHERE canon=?", (cn,))
        else:
            pct = int((nearest or {}).get("percent", 0) or 0)
            c.execute("INSERT INTO molecules(canon,smiles,first_player,first_target,first_seen,times_docked,props,synth,qed,drug_score,"
                      "drug_grade,nearest_name,nearest_pct,known_drug) VALUES(?,?,?,?,?,?,?,?,?,?,?,?,?,?)",
                      (cn, smiles, player, target_id, _now(), 1, json.dumps(props or {}), json.dumps(synth or {}), qed,
                       (prof or {}).get("score"), (prof or {}).get("grade_label"), (nearest or {}).get("name"), pct, 1 if pct >= 99 else 0))
        c.commit()
        return int(did)


def set_molecule_parent(canon_: str, parent_canon: str, parent_smiles: str) -> None:
    with _LOCK, closing(_conn()) as c:
        c.execute("UPDATE molecules SET parent_canon=?, parent_smiles=? WHERE canon=? AND (parent_canon IS NULL OR parent_canon='')",
                  (parent_canon, parent_smiles, canon_))
        c.commit()


def latest_docking(player: str, target_id: str, canon_: str, engine: str) -> Optional[Dict[str, Any]]:
    with _LOCK, closing(_conn()) as c:
        return _row(c.execute("SELECT * FROM dockings WHERE player=? AND target_id=? AND canon=? AND engine=? ORDER BY id DESC LIMIT 1",
                              (player, target_id, canon_, engine)).fetchone())


def best_docking(target_id: str, canon_: str, engine: str) -> Optional[Dict[str, Any]]:
    with _LOCK, closing(_conn()) as c:
        return _row(c.execute("SELECT * FROM dockings WHERE target_id=? AND canon=? AND engine=? ORDER BY energy ASC LIMIT 1",
                              (target_id, canon_, engine)).fetchone())


def list_molecules(limit: int = 100) -> List[Dict[str, Any]]:
    with _LOCK, closing(_conn()) as c:
        out = []
        for r in c.execute("SELECT * FROM molecules ORDER BY first_seen DESC LIMIT ?", (int(limit),)).fetchall():
            d = dict(r)
            for k in ("props", "synth"):
                try:
                    d[k] = json.loads(d[k] or "{}")
                except Exception:
                    d[k] = {}
            out.append(d)
        return out


# ── reference poses ──────────────────────────────────────────────────────────
def set_reference_pose(target_id, canon_, engine, energy, pose_pdb, source="blind") -> None:
    with _LOCK, closing(_conn()) as c:
        old = c.execute("SELECT energy FROM poses WHERE target_id=? AND canon=? AND engine=?", (target_id, canon_, engine)).fetchone()
        if old is None or float(energy) < old["energy"]:
            c.execute("INSERT OR REPLACE INTO poses VALUES(?,?,?,?,?,?,?)", (target_id, canon_, engine, float(energy), pose_pdb or "", source, _now()))
            c.commit()


def get_reference_pose(target_id, canon_, engine) -> Optional[Dict[str, Any]]:
    with _LOCK, closing(_conn()) as c:
        return _row(c.execute("SELECT * FROM poses WHERE target_id=? AND canon=? AND engine=?", (target_id, canon_, engine)).fetchone())


# ── records (the chain of "score to beat") ───────────────────────────────────
def head_record(target_id, engine, parent_canon) -> Optional[Dict[str, Any]]:
    with _LOCK, closing(_conn()) as c:
        return _row(c.execute("SELECT * FROM records WHERE target_id=? AND engine=? AND parent_canon=? ORDER BY id DESC LIMIT 1",
                              (target_id, engine, parent_canon)).fetchone())


def add_record(target_id, engine, parent_canon, holder_canon, holder_smiles, holder_name, player, energy, prev_energy, delta,
               rmsd_to_prev, decided_by, pose_pdb, docking_id) -> int:
    with _LOCK, closing(_conn()) as c:
        cur = c.execute("INSERT INTO records(target_id,engine,parent_canon,holder_canon,holder_smiles,holder_name,player,energy,"
                        "prev_energy,delta,rmsd_to_prev,decided_by,pose_pdb,docking_id,created) VALUES(?,?,?,?,?,?,?,?,?,?,?,?,?,?,?)",
                        (target_id, engine, parent_canon, holder_canon, holder_smiles, holder_name, player, float(energy),
                         None if prev_energy is None else float(prev_energy), delta, rmsd_to_prev, decided_by, pose_pdb or "",
                         docking_id, _now()))
        c.commit()
        return int(cur.lastrowid)


def record_history(target_id, engine, parent_canon="", limit=10) -> List[Dict[str, Any]]:
    q = "SELECT id,target_id,engine,parent_canon,holder_canon,holder_smiles,holder_name,player,energy,prev_energy,delta,rmsd_to_prev," \
        "decided_by,created FROM records WHERE target_id=? AND engine=?"
    args: list = [target_id, engine]
    if parent_canon:
        q += " AND parent_canon=?"; args.append(parent_canon)
    with _LOCK, closing(_conn()) as c:
        return [dict(r) for r in c.execute(q + " ORDER BY id DESC LIMIT ?", args + [int(limit)]).fetchall()]


# ── claims ───────────────────────────────────────────────────────────────────
def get_claim(docking_id: int) -> Optional[Dict[str, Any]]:
    with _LOCK, closing(_conn()) as c:
        r = c.execute("SELECT verdict FROM claims WHERE docking_id=?", (int(docking_id),)).fetchone()
        return json.loads(r["verdict"]) if r else None


def save_claim(docking_id: int, verdict: Dict[str, Any]) -> None:
    with _LOCK, closing(_conn()) as c:
        c.execute("INSERT OR REPLACE INTO claims VALUES(?,?,?)", (int(docking_id), json.dumps(verdict), _now()))
        c.commit()


# ── native (experimental) ligand pose ────────────────────────────────────────
def get_native(target_id: str) -> Optional[Dict[str, Any]]:
    with _LOCK, closing(_conn()) as c:
        return _row(c.execute("SELECT * FROM natives WHERE target_id=?", (target_id,)).fetchone())


def save_native(target_id: str, comp: str, pdb_lines: str) -> None:
    with _LOCK, closing(_conn()) as c:
        c.execute("INSERT OR REPLACE INTO natives VALUES(?,?,?,?)", (target_id, comp, pdb_lines, _now()))
        c.commit()
