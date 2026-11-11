"""Viewer-only endpoints for ligand previews and protonated receptor structures."""
import asyncio
import logging
from pathlib import Path

import numpy as np
from fastapi import APIRouter, HTTPException
from fastapi.responses import FileResponse
from pydantic import BaseModel
from rdkit import Chem

from .. import docking as dock_engine
from ..chem import smiles_to_mol3d
from ..config import RECEPTORS_DIR
from ..receptor_prep import RECEPTORS

logger = logging.getLogger("mobidock.viewer")
router = APIRouter()


class PreviewRequest(BaseModel):
    target_id: str
    smiles: str


def _preview_molecule(target_id: str, smiles: str) -> dict:
    box = dock_engine.BOX_PARAMS.get(target_id)
    if not box:
        raise ValueError(f"No docking box is available for target {target_id}.")
    mol = smiles_to_mol3d((smiles or "").strip())
    if mol is None or mol.GetNumConformers() == 0:
        raise ValueError("Could not build a 3-D conformer for this ligand.")

    center = np.asarray(box["center"], dtype=float)
    size = box["size"]
    conformer = mol.GetConformer()
    coordinates = conformer.GetPositions()
    shift = center - coordinates.mean(axis=0)
    shifted = coordinates + shift
    for index, point in enumerate(shifted):
        conformer.SetAtomPosition(index, point.tolist())

    heavy = shifted[[atom.GetAtomicNum() > 1 for atom in mol.GetAtoms()]]
    length = float(np.linalg.norm(heavy[:, None] - heavy[None], axis=2).max()) if len(heavy) > 1 else 0.0
    return {
        "sdf": Chem.MolToMolBlock(mol),
        "heavy_atoms": int(len(heavy)),
        "length": round(length, 1),
        "box_min_edge": round(float(min(size)), 1),
        "fits": length <= float(min(size)),
    }


@router.post("/api/insight/preview")
async def preview_ligand(req: PreviewRequest):
    """Return a centered conformer for the pre-docking protein viewer."""
    try:
        return await asyncio.to_thread(_preview_molecule, req.target_id, req.smiles)
    except ValueError as exc:
        raise HTTPException(status_code=422, detail=str(exc)) from exc
    except Exception as exc:
        logger.exception("Could not build ligand preview")
        raise HTTPException(status_code=422, detail="Could not build a 3-D model of this molecule") from exc


@router.get("/api/receptor/{target_id}/pdb_h")
async def protonated_receptor(target_id: str):
    """Serve the optional hydrogen-added receptor used by the 3-D viewer."""
    info = RECEPTORS.get(target_id)
    if not info:
        raise HTTPException(status_code=404, detail=f"Unknown target: {target_id}")
    clean = Path(info["pdb_path"]) if "pdb_path" in info else RECEPTORS_DIR / f"{target_id}_{info['pdb_id']}_clean.pdb"
    hydrogenated = clean.with_name(clean.name.replace("_clean.pdb", "_H.pdb"))
    if not hydrogenated.exists():
        raise HTTPException(status_code=404, detail="Protonated receptor file not found")
    return FileResponse(str(hydrogenated), media_type="text/plain")
