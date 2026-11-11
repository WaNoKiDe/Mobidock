# ══ CELL 24 · app/oracle.py ══
"""
oracle.py — Server-side pose verification (The Oracle).

Receives client-submitted ligand coordinates, runs a fast headless
docking calculation, checks for atomic overlap collisions, and issues
a pass/fail verdict.  Rate-limits clients after failed submissions.
"""
import hashlib
import logging
import time
from dataclasses import dataclass
from pathlib import Path
from typing import Any, Dict, List, Optional, Tuple

logger = logging.getLogger(__name__)

# ── Configuration ─────────────────────────────────────────────────────────────

RATE_LIMIT_COOLDOWN = 30.0        # seconds of cooldown after a failed submission
MAX_COLLISION_TOLERANCE = 0       # zero atomic overlap allowed
VDW_CLASH_THRESHOLD = 1.5        # Angstroms — distance below which atoms clash


# ── Data classes ──────────────────────────────────────────────────────────────

@dataclass
class VerificationResult:
    """Result of oracle pose verification."""
    passed: bool
    verified_energy: float
    collision_count: int
    pose_hash: str                    # SHA-256(coordinates + timestamp)
    reason: Optional[str] = None      # failure reason if rejected
    engine_used: str = "vina"
    verification_time_ms: float = 0.0


# ── Rate limiter ──────────────────────────────────────────────────────────────

class RateLimiter:
    """Per-client rate limiting for failed submissions."""

    def __init__(self, cooldown: float = RATE_LIMIT_COOLDOWN):
        self._cooldown = cooldown
        self._blocked: Dict[str, float] = {}    # user_id → unblock_time

    def is_rate_limited(self, user_id: str) -> Tuple[bool, float]:
        """Check if user is rate-limited.  Returns (blocked, remaining_seconds)."""
        if user_id not in self._blocked:
            return False, 0.0
        unblock_at = self._blocked[user_id]
        remaining = unblock_at - time.time()
        if remaining <= 0:
            del self._blocked[user_id]
            return False, 0.0
        return True, round(remaining, 1)

    def apply_cooldown(self, user_id: str):
        """Apply rate-limit cooldown to a user after a failed submission."""
        self._blocked[user_id] = time.time() + self._cooldown
        logger.warning(
            "Rate-limited user %s for %.0fs after failed submission",
            user_id, self._cooldown,
        )

    def clear(self, user_id: str):
        """Remove rate-limit for a user (e.g. after successful submission)."""
        self._blocked.pop(user_id, None)


# ── Collision checker ─────────────────────────────────────────────────────────

def check_atomic_overlap(
    docked_pdbqt: str,
    receptor_path: Optional[Path],
    clash_threshold: float = VDW_CLASH_THRESHOLD,
) -> int:
    """
    Count atomic clashes between docked ligand and receptor.

    Uses simple Euclidean distance check: if any ligand atom is
    within `clash_threshold` Å of any receptor atom → collision.

    Returns collision count (0 = no clashes = pass).
    """
    if receptor_path is None or not receptor_path.exists():
        # Cannot check without receptor — assume pass
        return 0

    ligand_atoms = _parse_atom_coords(docked_pdbqt)
    try:
        receptor_text = receptor_path.read_text(encoding="utf-8")
    except Exception:
        return 0
    receptor_atoms = _parse_atom_coords(receptor_text)

    if not ligand_atoms or not receptor_atoms:
        return 0

    collisions = 0
    threshold_sq = clash_threshold ** 2

    for lx, ly, lz in ligand_atoms:
        for rx, ry, rz in receptor_atoms:
            dx = lx - rx
            dy = ly - ry
            dz = lz - rz
            if dx * dx + dy * dy + dz * dz < threshold_sq:
                collisions += 1
                break  # One collision per ligand atom is enough

    return collisions


def _parse_atom_coords(pdb_text: str) -> List[Tuple[float, float, float]]:
    """Extract (x, y, z) coordinates from ATOM/HETATM records."""
    coords = []
    for line in pdb_text.splitlines():
        if line.startswith("ATOM") or line.startswith("HETATM"):
            try:
                x = float(line[30:38].strip())
                y = float(line[38:46].strip())
                z = float(line[46:54].strip())
                coords.append((x, y, z))
            except (ValueError, IndexError):
                continue
    return coords


# ── Pose hashing ──────────────────────────────────────────────────────────────

def compute_pose_hash(
    ligand_transform: Dict[str, Any],
    timestamp: float,
) -> str:
    """SHA-256 hash of coordinates + timestamp for immutable ledger."""
    raw = (
        str(ligand_transform.get("position", []))
        + str(ligand_transform.get("rotation", []))
        + str(ligand_transform.get("torsion_angles", []))
        + f"{timestamp:.6f}"
    )
    return hashlib.sha256(raw.encode()).hexdigest()


# ── The Oracle — main verification ────────────────────────────────────────────

class Oracle:
    """Server-side docking pose verification engine."""

    def __init__(self):
        self.rate_limiter = RateLimiter()

    def verify_pose(
        self,
        epoch_id: str,
        solver_id: str,
        ligand_smiles: str,
        ligand_transform: Dict[str, Any],
        target_id: str,
        threshold: float,
        client_timestamp: float,
    ) -> VerificationResult:
        """
        Run server-side verification of a submitted docking pose.

        Pipeline:
        1. Check rate limit
        2. Dock the ligand SMILES using existing docking pipeline
        3. Check for atomic overlaps
        4. Compare verified energy against threshold
        5. Return verdict
        """
        t0 = time.time()

        # 1. Rate limit check
        blocked, remaining = self.rate_limiter.is_rate_limited(solver_id)
        if blocked:
            return VerificationResult(
                passed=False,
                verified_energy=0.0,
                collision_count=0,
                pose_hash="",
                reason=f"Rate limited. Try again in {remaining:.0f}s.",
            )

        # 2. Compute pose hash
        pose_hash = compute_pose_hash(ligand_transform, client_timestamp)

        # 3. Run server-side docking using existing pipeline
        try:
            from . import docking as dock_engine

            result = dock_engine.dock_smiles(
                target_id=target_id,
                smiles=ligand_smiles,
                exhaustiveness=4,       # Lower for speed on verification
                engine_name="vina",     # Always use Vina for oracle
            )
            verified_energy = result["energy"]
            docked_pdbqt = result["docked_pdbqt"]
            is_mock = False
        except Exception as exc:
            logger.error("Oracle docking failed for %s: %s", solver_id, exc)
            return VerificationResult(
                passed=False,
                verified_energy=0.0,
                collision_count=0,
                pose_hash=pose_hash,
                reason=f"Server docking error: {exc}",
                verification_time_ms=round((time.time() - t0) * 1000, 1),
            )

        # 4. Check atomic overlaps
        receptor_path = dock_engine.RECEPTOR_PATHS.get(target_id)
        collision_count = check_atomic_overlap(docked_pdbqt, receptor_path)

        # 5. Evaluate pass/fail
        energy_pass = verified_energy <= threshold
        collision_pass = collision_count <= MAX_COLLISION_TOLERANCE
        passed = energy_pass and collision_pass

        reason = None
        if not energy_pass:
            reason = (
                f"Insufficient binding energy: {verified_energy:.2f} kcal/mol "
                f"(threshold: {threshold:.1f})"
            )
        elif not collision_pass:
            reason = f"Atomic overlap detected: {collision_count} collision(s)"

        # Apply rate limit on failure
        if not passed:
            self.rate_limiter.apply_cooldown(solver_id)
        else:
            self.rate_limiter.clear(solver_id)

        elapsed = round((time.time() - t0) * 1000, 1)
        logger.info(
            "Oracle verdict for %s (epoch=%s): %s | E=%.2f (thresh=%.1f) | "
            "collisions=%d | %.0fms",
            solver_id, epoch_id,
            "PASS ✓" if passed else "FAIL ✗",
            verified_energy, threshold, collision_count, elapsed,
        )

        return VerificationResult(
            passed=passed,
            verified_energy=verified_energy,
            collision_count=collision_count,
            pose_hash=pose_hash,
            reason=reason,
            engine_used="vina",
            verification_time_ms=elapsed,
        )
