# ══ CELL 23 · app/epoch.py ══
"""
epoch.py — Epoch-based target distribution system.

Manages the lifecycle of "Target Epochs": timed competitive docking
challenges broadcast to all connected players.  Each epoch selects a
protein target + reference ligand, sets a difficulty threshold, and
offers a bounty of Compute Credits.
"""
import asyncio
import hashlib
import json
import logging
import random
import time
from dataclasses import dataclass, field, asdict
from pathlib import Path
from typing import Any, Dict, List, Optional, Callable

logger = logging.getLogger(__name__)

# ── Epoch configuration ──────────────────────────────────────────────────────

DEFAULT_EPOCH_DURATION = 300        # 5 minutes per epoch
BASE_BOUNTY            = 500        # Compute Credits for first-to-solve
MIN_THRESHOLD          = -12.0      # kcal/mol — hardest difficulty
MAX_THRESHOLD          = -5.0       # kcal/mol — easiest difficulty

# Pre-seeded ligand pool per target family
LIGAND_POOLS: Dict[str, List[Dict[str, str]]] = {
    "general": [
        {"name": "Aspirin",        "smiles": "CC(=O)Oc1ccccc1C(=O)O"},
        {"name": "Ibuprofen",      "smiles": "CC(C)Cc1ccc(cc1)C(C)C(=O)O"},
        {"name": "Caffeine",       "smiles": "Cn1c(=O)c2c(ncn2C)n(c1=O)C"},
        {"name": "Dopamine",       "smiles": "NCCc1ccc(O)c(O)c1"},
        {"name": "Serotonin",      "smiles": "NCCc1c[nH]c2ccc(O)cc12"},
        {"name": "Melatonin",      "smiles": "COc1ccc2[nH]cc(CCNC(C)=O)c2c1"},
        {"name": "Resveratrol",    "smiles": "Oc1ccc(cc1)/C=C/c1cc(O)cc(O)c1"},
        {"name": "Paracetamol",    "smiles": "CC(=O)Nc1ccc(O)cc1"},
        {"name": "Naproxen",       "smiles": "COc1ccc2cc([C@H](C)C(=O)O)ccc2c1"},
        {"name": "Donepezil",      "smiles": "COC1=C(C=C2C(=C1)CC(C2=O)CC3CCN(CC3)CC4=CC=CC=C4)OC"},
        {"name": "Selegiline",     "smiles": "C#CCN(C)[C@H](C)Cc1ccccc1"},
        {"name": "Quercetin",      "smiles": "Oc1cc(O)c2c(c1)oc(-c1ccc(O)c(O)c1)c(O)c2=O"},
    ],
}


# ── Data classes ──────────────────────────────────────────────────────────────

@dataclass
class TargetEpoch:
    """A single competitive docking epoch."""
    epoch_id: str
    protein_id: str                    # target_id from RECEPTORS
    ligand_smiles: str
    ligand_name: str
    difficulty_threshold_kcal_mol: float
    bounty_amount: int
    broadcast_timestamp: float
    epoch_timer: int                   # seconds
    start_time: float = field(default_factory=time.time)
    solved_by: Optional[str] = None
    solved_at: Optional[float] = None

    def remaining_seconds(self) -> float:
        return max(0.0, self.epoch_timer - (time.time() - self.start_time))

    def is_expired(self) -> bool:
        return self.remaining_seconds() <= 0

    def to_broadcast(self) -> Dict[str, Any]:
        """Payload sent to WebSocket clients."""
        return {
            "epoch_id": self.epoch_id,
            "protein_id": self.protein_id,
            "ligand_smiles": self.ligand_smiles,
            "ligand_name": self.ligand_name,
            "difficulty_threshold_kcal_mol": self.difficulty_threshold_kcal_mol,
            "bounty_amount": self.bounty_amount,
            "broadcast_timestamp": self.broadcast_timestamp,
            "epoch_timer": self.epoch_timer,
            "remaining_seconds": round(self.remaining_seconds(), 1),
            "solved_by": self.solved_by,
            "solved_at": self.solved_at,
        }


# ── Epoch Manager ─────────────────────────────────────────────────────────────

class EpochManager:
    """
    Manages the active epoch and handles rotation.

    Usage (inside FastAPI lifespan):
        manager = EpochManager(target_ids=["BACE1", "MAOA", ...])
        asyncio.create_task(manager.run())
    """

    def __init__(
        self,
        target_ids: List[str],
        epoch_duration: int = DEFAULT_EPOCH_DURATION,
        base_bounty: int = BASE_BOUNTY,
    ):
        self.target_ids = list(target_ids) if target_ids else ["BACE1"]
        self.epoch_duration = epoch_duration
        self.base_bounty = base_bounty
        self._current: Optional[TargetEpoch] = None
        self._epoch_counter = 0
        self._history: List[Dict] = []
        self._subscribers: List[Callable] = []  # WebSocket broadcast callbacks
        self._running = False
        self._rng = random.Random(42)

    # ── Public API ────────────────────────────────────────────────────────

    @property
    def current_epoch(self) -> Optional[TargetEpoch]:
        return self._current

    def get_active_epoch(self) -> Optional[Dict[str, Any]]:
        """Return the active epoch as a dict, or None."""
        if self._current is None:
            return None
        return self._current.to_broadcast()

    def is_epoch_active(self, epoch_id: str) -> bool:
        """Check if a specific epoch is the current active one."""
        if self._current is None:
            return False
        return (
            self._current.epoch_id == epoch_id
            and not self._current.is_expired()
        )

    def mark_solved(self, epoch_id: str, solver_id: str) -> bool:
        """Mark the current epoch as solved by a player.  Returns True if this
        is the *first* solve (epoch winner)."""
        if not self.is_epoch_active(epoch_id):
            return False
        if self._current.solved_by is None:
            self._current.solved_by = solver_id
            self._current.solved_at = time.time()
            logger.info(
                "Epoch %s FIRST SOLVE by %s (%.1fs into epoch)",
                epoch_id, solver_id,
                time.time() - self._current.start_time,
            )
            return True
        return False  # Already solved — subsequent submission

    def subscribe(self, callback: Callable):
        """Register a callback for epoch broadcasts."""
        self._subscribers.append(callback)

    def unsubscribe(self, callback: Callable):
        """Unregister a broadcast callback."""
        if callback in self._subscribers:
            self._subscribers.remove(callback)

    # ── Epoch lifecycle ───────────────────────────────────────────────────

    def advance_epoch(self) -> TargetEpoch:
        """Create and activate a new epoch."""
        # Archive current
        if self._current is not None:
            self._history.append(self._current.to_broadcast())

        self._epoch_counter += 1
        epoch_id = f"EP-{self._epoch_counter:04d}"

        # Round-robin through targets with some randomness
        target_idx = (self._epoch_counter - 1) % len(self.target_ids)
        protein_id = self.target_ids[target_idx]

        # Pick a random ligand from the pool
        pool = LIGAND_POOLS.get(protein_id, LIGAND_POOLS["general"])
        ligand = self._rng.choice(pool)

        # Scale difficulty based on epoch number (gets harder over time)
        progress = min(1.0, self._epoch_counter / 20.0)
        threshold = round(
            MAX_THRESHOLD + (MIN_THRESHOLD - MAX_THRESHOLD) * progress * 0.5,
            1,
        )

        # Bounty increases slightly with difficulty
        bounty = self.base_bounty + int(50 * progress * self._epoch_counter % 10)

        epoch = TargetEpoch(
            epoch_id=epoch_id,
            protein_id=protein_id,
            ligand_smiles=ligand["smiles"],
            ligand_name=ligand["name"],
            difficulty_threshold_kcal_mol=threshold,
            bounty_amount=bounty,
            broadcast_timestamp=time.time(),
            epoch_timer=self.epoch_duration,
        )
        self._current = epoch
        logger.info(
            "New epoch %s: target=%s, ligand=%s, threshold=%.1f, bounty=%d",
            epoch_id, protein_id, ligand["name"], threshold, bounty,
        )
        return epoch

    async def _broadcast(self, event_type: str, data: Dict):
        """Send event to all WebSocket subscribers."""
        message = {"type": event_type, "data": data}
        dead = []
        for cb in self._subscribers:
            try:
                await cb(message)
            except Exception:
                dead.append(cb)
        for cb in dead:
            self._subscribers.remove(cb)

    async def run(self):
        """Main loop — advance epochs and broadcast state."""
        self._running = True
        logger.info("EpochManager started (duration=%ds)", self.epoch_duration)

        while self._running:
            # Advance to new epoch
            epoch = self.advance_epoch()
            await self._broadcast("new_epoch", epoch.to_broadcast())

            # Tick every second until epoch expires
            while not epoch.is_expired() and self._running:
                await asyncio.sleep(1.0)

                # Broadcast periodic timer updates (every 5 seconds)
                remaining = epoch.remaining_seconds()
                if int(remaining) % 5 == 0 and remaining > 0:
                    await self._broadcast("epoch_tick", {
                        "epoch_id": epoch.epoch_id,
                        "remaining_seconds": round(remaining, 1),
                        "solved_by": epoch.solved_by,
                    })

                # If solved, broadcast the solve event once
                if epoch.solved_by and epoch.solved_at:
                    # Only broadcast once — check if we already did
                    if not hasattr(epoch, '_solve_broadcast'):
                        epoch._solve_broadcast = True
                        await self._broadcast("epoch_solved", {
                            "epoch_id": epoch.epoch_id,
                            "solved_by": epoch.solved_by,
                            "solved_at": epoch.solved_at,
                            "remaining_seconds": round(remaining, 1),
                        })

            # Epoch ended
            await self._broadcast("epoch_expired", {
                "epoch_id": epoch.epoch_id,
                "solved_by": epoch.solved_by,
            })

        logger.info("EpochManager stopped.")

    def stop(self):
        """Signal the run loop to exit."""
        self._running = False
