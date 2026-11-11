# ══ CELL 07 · app/runtime.py ══
"""
runtime.py — Shared singletons and request/response models used by several route files.

main.py fills these in during startup (lifespan). Route files do:
    from .. import runtime as rt
and then use rt.epoch_manager, rt.oracle, rt.entropy_engine, rt.ws_connections.
"""
from typing import List, Optional

from pydantic import BaseModel, Field

epoch_manager = None      # EpochManager   (set in main.lifespan)
oracle = None             # Oracle
entropy_engine = None     # EntropyEngine
ws_connections: List = [] # connected WebSocket clients


# ── Docking models ────────────────────────────────────────────────────────────
class DockRequest(BaseModel):
    target_id: str
    smiles: str
    player_name: Optional[str] = "Anonymous"
    exhaustiveness: Optional[int] = 8
    n_poses: Optional[int] = 10
    engine: Optional[str] = "vina"
    lipinski_mode: str = "standard"
    synth_mode: str = "standard"
    box: Optional[dict] = None
    diffdock_api_key: Optional[str] = None


class DockCompareRequest(BaseModel):
    target_id: str
    smiles: str
    player_name: Optional[str] = "Anonymous"
    engine: Optional[str] = None  # None runs all registered engines; a name runs one engine
    exhaustiveness: Optional[int] = 8
    n_poses: Optional[int] = 10
    lipinski_mode: str = "standard"
    synth_mode: str = "standard"
    ligand_source: str = "chosen"
    ligand_label: str = "Chosen ligand"
    box: Optional[dict] = None
    diffdock_api_key: Optional[str] = None


class DockResponse(BaseModel):
    target_id: str
    smiles: str
    player_name: str
    energy: float
    qed: float
    score: float
    rank: int
    props: dict
    docked_pdbqt: str
    docked_pdb: str = ""          # best pose as PDB — what the 3D viewer renders
    ligand_pdbqt: str
    is_mock: bool = False
    engine: str
    engine_metadata: dict
    box: dict = Field(default_factory=dict)  # grid box used for this run
    synth: dict = Field(default_factory=dict)   # synthesizability report


class BoxUpdateRequest(BaseModel):
    center: List[float]           # [x, y, z]
    size: List[float]             # [sx, sy, sz]  full edge lengths of the grid box (Å)


class PdbImportRequest(BaseModel):
    pdb_id: str


class EngineKeyRequest(BaseModel):
    api_key: str


# ── Game models (epoch / ledger / entropy) ────────────────────────────────────
class LigandTransform(BaseModel):
    position: List[float] = Field(default_factory=lambda: [0.0, 0.0, 0.0])
    rotation: List[float] = Field(default_factory=lambda: [0.0, 0.0, 0.0, 1.0])
    torsion_angles: List[float] = Field(default_factory=list)


class SolverSubmission(BaseModel):
    epoch_id: str
    solver_id: str = "Anonymous"
    client_timestamp: float = 0.0
    ligand_transform: LigandTransform = Field(default_factory=LigandTransform)


class EntropySubmitRequest(BaseModel):
    user_id: str
    epoch_id: str
    target_id: str
    pose_hash: str
    ligand_smiles: str
    binding_energy: float
