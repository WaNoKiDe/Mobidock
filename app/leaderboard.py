# ══ CELL 21 · app/leaderboard.py ══
"""
leaderboard.py  —  Thread-safe JSON-backed leaderboard for MedChem Game
"""
import json
import threading
from pathlib import Path
from typing import List, Dict, Any, Optional

from .synth import composite_score

LEADERBOARD_FILE = Path(__file__).parent.parent / "data" / "leaderboard.json"
_lock = threading.Lock()

# ── Initialise empty file if needed ──────────────────────────────────────────

def _load() -> Dict[str, List]:
    if not LEADERBOARD_FILE.exists():
        LEADERBOARD_FILE.parent.mkdir(parents=True, exist_ok=True)
        return {}
    try:
        with open(LEADERBOARD_FILE, "r", encoding="utf-8") as f:
            return json.load(f)
    except (json.JSONDecodeError, OSError):
        return {}


def _save(data: Dict[str, List]) -> None:
    LEADERBOARD_FILE.parent.mkdir(parents=True, exist_ok=True)
    with open(LEADERBOARD_FILE, "w", encoding="utf-8") as f:
        json.dump(data, f, indent=2)


# ── Public API ────────────────────────────────────────────────────────────────

def add_score(target_id: str, player_name: str, energy: float, qed: float,
              synth_norm: Optional[float] = None, sa: Optional[float] = None) -> int:
    """Add a score entry and return the player's rank (1-indexed) for that target."""
    with _lock:
        data = _load()
        if target_id not in data:
            data[target_id] = []

        entry = {
            "player": player_name,
            "energy": round(energy, 3),
            "qed": round(qed, 3),
            "synth": sa,                      # SA score (1 easy .. 10 very hard); None = not scored
            # composite score: lower energy, higher QED and easier synthesis all score higher
            "score": composite_score(energy, qed, synth_norm),
        }
        data[target_id].append(entry)

        # Sort: higher composite score = better
        data[target_id].sort(key=lambda x: x["score"], reverse=True)

        # Keep top 100 per target
        data[target_id] = data[target_id][:100]

        _save(data)

        # Find rank
        for idx, e in enumerate(data[target_id]):
            if e["player"] == player_name and e["energy"] == entry["energy"]:
                return idx + 1
        return len(data[target_id])


def get_leaderboard(target_id: str = None, top_n: int = 10) -> Dict[str, Any]:
    """Return top_n entries per target (or all targets if target_id is None)."""
    with _lock:
        data = _load()

    if target_id:
        entries = data.get(target_id, [])[:top_n]
        return {target_id: entries}

    return {tid: entries[:top_n] for tid, entries in data.items()}
