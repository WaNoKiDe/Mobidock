"""Route-specific RDKit descriptor warnings for the loaded ligand."""
import asyncio
import logging

from fastapi import APIRouter, HTTPException
from pydantic import BaseModel
from rdkit import Chem

from ..drug_warnings import warnings_all

logger = logging.getLogger("mobidock.drug_warnings")
router = APIRouter()


class DrugWarningsRequest(BaseModel):
    smiles: str


def _warnings(smiles: str) -> dict:
    mol = Chem.MolFromSmiles((smiles or "").strip())
    if mol is None:
        raise ValueError("That SMILES is not valid")
    return warnings_all(mol)


@router.post("/api/drug-warnings")
async def drug_warnings(req: DrugWarningsRequest):
    try:
        return await asyncio.to_thread(_warnings, req.smiles)
    except ValueError as exc:
        raise HTTPException(status_code=422, detail=str(exc))
    except Exception:
        logger.exception("Drug warning calculation failed")
        raise HTTPException(status_code=500, detail="Could not calculate RDKit property warnings")
