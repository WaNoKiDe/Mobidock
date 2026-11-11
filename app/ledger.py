# ══ CELL 22 · app/ledger.py ══
"""
ledger.py — Append-only resolution ledger & Compute Credits economy.

Every successful epoch solve is recorded as an immutable entry.
Credits are derived from the ledger — the JSON credits file is a
materialised view rebuilt on startup.
"""
import json
import math
import logging
import threading
import time
from dataclasses import dataclass, asdict
from pathlib import Path
from typing import Any, Dict, List, Optional

logger = logging.getLogger(__name__)

# ── Paths ─────────────────────────────────────────────────────────────────────

DATA_DIR = Path(__file__).parent.parent / "data"
LEDGER_FILE = DATA_DIR / "ledger.jsonl"
CREDITS_FILE = DATA_DIR / "credits.json"

_lock = threading.Lock()


# ── Data classes ──────────────────────────────────────────────────────────────

@dataclass
class ResolutionEntry:
    """A single immutable ledger entry."""
    epoch_id: str
    solver_user_id: str
    docked_pose_hash: str
    verified_binding_energy: float
    credits_awarded: int
    timestamp: float
    is_first_solve: bool
    target_id: str
    ligand_smiles: str
    entry_type: str = "epoch_solve"    # epoch_solve | entropy_debit | entropy_bonus


# ── Reward calculation ────────────────────────────────────────────────────────

def calculate_reward(
    bounty: int,
    is_first_solve: bool,
    time_delta_seconds: float,
) -> int:
    """
    Calculate credit reward for a solve.

    - First-to-solve: 100% of bounty
    - Subsequent: time-decayed  bounty * exp(-0.05 * delta)
    - Minimum: 10 credits (always worth submitting)
    """
    if is_first_solve:
        return bounty

    decayed = bounty * math.exp(-0.05 * time_delta_seconds)
    return max(10, int(round(decayed)))


# ── Ledger I/O ────────────────────────────────────────────────────────────────

def _ensure_data_dir():
    DATA_DIR.mkdir(parents=True, exist_ok=True)


def append_entry(entry: ResolutionEntry) -> None:
    """Append a single entry to the immutable JSONL ledger."""
    _ensure_data_dir()
    with _lock:
        with open(LEDGER_FILE, "a", encoding="utf-8") as f:
            f.write(json.dumps(asdict(entry)) + "\n")
        # Update materialised credit view
        _update_credits(entry.solver_user_id, entry.credits_awarded)

    logger.info(
        "Ledger entry: %s solved %s (%+d credits, first=%s)",
        entry.solver_user_id, entry.epoch_id,
        entry.credits_awarded, entry.is_first_solve,
    )


def load_all_entries() -> List[Dict[str, Any]]:
    """Read the full ledger (for startup rebuild or admin)."""
    if not LEDGER_FILE.exists():
        return []
    entries = []
    with open(LEDGER_FILE, "r", encoding="utf-8") as f:
        for line in f:
            line = line.strip()
            if line:
                try:
                    entries.append(json.loads(line))
                except json.JSONDecodeError:
                    continue
    return entries


# ── Credits ───────────────────────────────────────────────────────────────────

def _load_credits() -> Dict[str, int]:
    """Load materialised credit balances."""
    if not CREDITS_FILE.exists():
        return {}
    try:
        with open(CREDITS_FILE, "r", encoding="utf-8") as f:
            return json.load(f)
    except (json.JSONDecodeError, OSError):
        return {}


def _save_credits(credits: Dict[str, int]) -> None:
    _ensure_data_dir()
    with open(CREDITS_FILE, "w", encoding="utf-8") as f:
        json.dump(credits, f, indent=2)


def _update_credits(user_id: str, delta: int) -> None:
    """Atomically update a user's credit balance."""
    credits = _load_credits()
    credits[user_id] = credits.get(user_id, 0) + delta
    _save_credits(credits)


def get_player_credits(user_id: str) -> int:
    """Return current credit balance for a player."""
    with _lock:
        credits = _load_credits()
        return credits.get(user_id, 0)


def deduct_credits(user_id: str, amount: int) -> bool:
    """
    Deduct credits from a player.  Returns True if successful,
    False if insufficient balance.
    """
    with _lock:
        credits = _load_credits()
        current = credits.get(user_id, 0)
        if current < amount:
            return False
        credits[user_id] = current - amount
        _save_credits(credits)
        return True


def get_player_history(user_id: str, limit: int = 50) -> List[Dict[str, Any]]:
    """Return resolution entries for a specific player (newest first)."""
    all_entries = load_all_entries()
    player_entries = [
        e for e in all_entries if e.get("solver_user_id") == user_id
    ]
    return sorted(
        player_entries,
        key=lambda e: e.get("timestamp", 0),
        reverse=True,
    )[:limit]


def get_all_credit_balances() -> Dict[str, int]:
    """Return all player credit balances (for dashboard)."""
    with _lock:
        return _load_credits()


def rebuild_credits_from_ledger() -> Dict[str, int]:
    """Rebuild the materialised credit view from the full ledger."""
    entries = load_all_entries()
    credits: Dict[str, int] = {}
    for e in entries:
        uid = e.get("solver_user_id", "")
        delta = e.get("credits_awarded", 0)
        credits[uid] = credits.get(uid, 0) + delta

    with _lock:
        _save_credits(credits)

    logger.info(
        "Credits rebuilt from ledger: %d players, %d total entries",
        len(credits), len(entries),
    )
    return credits
