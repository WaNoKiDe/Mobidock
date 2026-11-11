# ══ CELL 25 · app/entropy_engine.py ══
"""
entropy_engine.py — Thermodynamic stress-testing (Entropy Engine).

Players spend Compute Credits to submit winning docking poses to an
asynchronous Molecular Dynamics (MD) simulation queue.  The engine
measures RMSD over simulated timesteps to classify poses as
"Thermodynamically Verified" (stable) or "False Equilibrium" (unstable).

NOTE: Real MD simulation requires OpenMM/GROMACS + GPU.  This module
ships with a physically-plausible mock that generates RMSD trajectories
based on binding energy and structural noise.  The interface is designed
so a real MD backend can be plugged in without changing API contracts.
"""
import asyncio
import hashlib
import json
import logging
import math
import random
import time
import uuid
from dataclasses import dataclass, field, asdict
from enum import Enum
from pathlib import Path
from typing import Any, Dict, List, Optional

logger = logging.getLogger(__name__)

# ── Configuration ─────────────────────────────────────────────────────────────

ENTROPY_COST_CREDITS = 100        # Credits to run one MD test
RMSD_STABLE_THRESHOLD = 2.0      # Angstroms — below this = stable
SCORE_MULTIPLIER = 3.0            # 3× multiplier for verified poses
MD_SIM_TIMESTEPS = 100            # Number of simulated MD steps
SIM_DELAY_SECONDS = (5, 15)       # (min, max) simulated processing time

DATA_DIR = Path(__file__).parent.parent / "data"
JOBS_FILE = DATA_DIR / "entropy_jobs.json"


# ── Enums / Data classes ──────────────────────────────────────────────────────

class JobStatus(str, Enum):
    QUEUED    = "queued"
    RUNNING   = "running"
    COMPLETE  = "complete"
    FAILED    = "failed"


class Verdict(str, Enum):
    VERIFIED  = "thermodynamically_verified"
    REJECTED  = "false_equilibrium"
    PENDING   = "pending"


@dataclass
class EntropyJob:
    """Represents a single MD stress-test job."""
    job_id: str
    user_id: str
    epoch_id: str
    target_id: str
    pose_hash: str
    ligand_smiles: str
    binding_energy: float
    credits_cost: int
    status: str = JobStatus.QUEUED
    verdict: str = Verdict.PENDING
    rmsd_trajectory: List[float] = field(default_factory=list)
    mean_rmsd: float = 0.0
    max_rmsd: float = 0.0
    multiplier_awarded: float = 0.0
    credits_bonus: int = 0
    created_at: float = field(default_factory=time.time)
    completed_at: Optional[float] = None
    error_message: Optional[str] = None


# ── Job persistence ───────────────────────────────────────────────────────────

def _load_jobs() -> Dict[str, Dict]:
    if not JOBS_FILE.exists():
        return {}
    try:
        with open(JOBS_FILE, "r", encoding="utf-8") as f:
            return json.load(f)
    except (json.JSONDecodeError, OSError):
        return {}


def _save_jobs(jobs: Dict[str, Dict]) -> None:
    DATA_DIR.mkdir(parents=True, exist_ok=True)
    with open(JOBS_FILE, "w", encoding="utf-8") as f:
        json.dump(jobs, f, indent=2)


def _persist_job(job: EntropyJob) -> None:
    jobs = _load_jobs()
    jobs[job.job_id] = asdict(job)
    _save_jobs(jobs)


# ── Mock MD simulation ────────────────────────────────────────────────────────

def _simulate_md_trajectory(
    binding_energy: float,
    timesteps: int = MD_SIM_TIMESTEPS,
    seed: Optional[int] = None,
) -> List[float]:
    """
    Generate a physically-plausible RMSD trajectory.

    Model:
    - Stronger binding (more negative energy) → lower baseline RMSD
    - Each timestep adds Brownian noise + thermal drift
    - Unstable poses show monotonic RMSD increase (dissociation)
    - Stable poses oscillate around a low baseline

    Parameters
    ----------
    binding_energy : float
        kcal/mol (negative = stronger binding)
    timesteps : int
        Number of simulation steps

    Returns
    -------
    List of RMSD values (Angstroms) at each timestep
    """
    rng = random.Random(seed or hash(str(binding_energy) + str(time.time())))

    # Baseline RMSD inversely proportional to binding strength
    # Strong binder (E = -10): baseline ~ 0.5 Å
    # Weak binder (E = -3):  baseline ~ 2.5 Å
    baseline = max(0.3, 0.5 + (binding_energy + 10.0) * 0.3)

    # Stability factor: probability of staying near baseline
    # More negative energy → more stable
    stability = min(1.0, max(0.0, (-binding_energy - 3.0) / 8.0))

    trajectory = []
    current_rmsd = 0.0

    for step in range(timesteps):
        t_fraction = step / max(1, timesteps - 1)

        # Thermal noise (Brownian motion)
        noise = rng.gauss(0, 0.15)

        # Drift: stable poses stay near baseline, unstable ones drift up
        if stability > 0.6:
            # Stable: oscillate around baseline
            drift = (baseline - current_rmsd) * 0.1
        else:
            # Unstable: gradual dissociation with acceleration
            drift = 0.02 + 0.03 * t_fraction + rng.uniform(0, 0.05)

        # Temperature fluctuation events (occasional spikes)
        if rng.random() < 0.05:
            noise += rng.uniform(0.2, 0.8) * (1.0 - stability)

        current_rmsd = max(0.0, current_rmsd + drift + noise)
        trajectory.append(round(current_rmsd, 3))

    return trajectory


# ── Entropy Engine ────────────────────────────────────────────────────────────

class EntropyEngine:
    """
    Manages the thermodynamic stress-test job queue.

    Usage:
        engine = EntropyEngine()
        job = engine.submit_job(user_id, epoch_id, ...)
        # Job runs asynchronously
        status = engine.get_job(job.job_id)
    """

    def __init__(self):
        self._active_jobs: Dict[str, EntropyJob] = {}
        self._load_persisted_jobs()

    def _load_persisted_jobs(self):
        """Reload completed jobs from disk."""
        jobs = _load_jobs()
        for jid, jdata in jobs.items():
            try:
                self._active_jobs[jid] = EntropyJob(**{
                    k: v for k, v in jdata.items()
                    if k in EntropyJob.__dataclass_fields__
                })
            except Exception:
                continue

    def submit_job(
        self,
        user_id: str,
        epoch_id: str,
        target_id: str,
        pose_hash: str,
        ligand_smiles: str,
        binding_energy: float,
        credits_cost: int = ENTROPY_COST_CREDITS,
    ) -> EntropyJob:
        """
        Create a new MD stress-test job.

        NOTE: Credit deduction must be done by the caller BEFORE
        calling this method (via ledger.deduct_credits).
        """
        job = EntropyJob(
            job_id=f"MD-{uuid.uuid4().hex[:8].upper()}",
            user_id=user_id,
            epoch_id=epoch_id,
            target_id=target_id,
            pose_hash=pose_hash,
            ligand_smiles=ligand_smiles,
            binding_energy=binding_energy,
            credits_cost=credits_cost,
        )
        self._active_jobs[job.job_id] = job
        _persist_job(job)

        logger.info(
            "Entropy job %s queued: user=%s, E=%.2f, cost=%d CC",
            job.job_id, user_id, binding_energy, credits_cost,
        )
        return job

    async def process_job(self, job_id: str) -> EntropyJob:
        """
        Execute the MD simulation for a queued job.
        Runs asynchronously with simulated delay.
        """
        job = self._active_jobs.get(job_id)
        if job is None:
            raise ValueError(f"Job {job_id} not found")

        job.status = JobStatus.RUNNING
        _persist_job(job)

        try:
            # Simulated processing delay
            delay = random.uniform(*SIM_DELAY_SECONDS)
            await asyncio.sleep(delay)

            # Run mock MD
            trajectory = _simulate_md_trajectory(
                binding_energy=job.binding_energy,
                timesteps=MD_SIM_TIMESTEPS,
            )
            job.rmsd_trajectory = trajectory
            job.mean_rmsd = round(sum(trajectory) / len(trajectory), 3)
            job.max_rmsd = round(max(trajectory), 3)

            # Verdict
            if job.mean_rmsd < RMSD_STABLE_THRESHOLD:
                job.verdict = Verdict.VERIFIED
                job.multiplier_awarded = SCORE_MULTIPLIER
                job.credits_bonus = int(job.credits_cost * (SCORE_MULTIPLIER - 1))
                logger.info(
                    "Entropy job %s: VERIFIED ✓ (mean RMSD=%.2f Å, bonus=%d CC)",
                    job.job_id, job.mean_rmsd, job.credits_bonus,
                )
            else:
                job.verdict = Verdict.REJECTED
                job.multiplier_awarded = 0.0
                job.credits_bonus = 0
                logger.info(
                    "Entropy job %s: REJECTED ✗ (mean RMSD=%.2f Å — false equilibrium)",
                    job.job_id, job.mean_rmsd,
                )

            job.status = JobStatus.COMPLETE
            job.completed_at = time.time()

        except Exception as exc:
            job.status = JobStatus.FAILED
            job.error_message = str(exc)
            job.completed_at = time.time()
            logger.error("Entropy job %s failed: %s", job.job_id, exc)

        _persist_job(job)
        return job

    def get_job(self, job_id: str) -> Optional[Dict[str, Any]]:
        """Get a job's current state."""
        job = self._active_jobs.get(job_id)
        if job is None:
            # Try loading from disk
            jobs = _load_jobs()
            return jobs.get(job_id)
        return asdict(job)

    def get_player_jobs(self, user_id: str, limit: int = 20) -> List[Dict[str, Any]]:
        """Get all jobs for a player (newest first)."""
        all_jobs = []
        for job in self._active_jobs.values():
            if job.user_id == user_id:
                all_jobs.append(asdict(job))

        return sorted(
            all_jobs,
            key=lambda j: j.get("created_at", 0),
            reverse=True,
        )[:limit]
