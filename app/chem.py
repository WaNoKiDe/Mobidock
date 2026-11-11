# ══ CELL 06 · app/chem.py ══
"""
chem.py — Ligand chemistry helpers (RDKit + Meeko). No docking logic here.

SMILES → 3D molecule → PDBQT / SDF, plus QED and property calculations.
"""
import logging
import math
from typing import Any, Dict, Optional

from rdkit import Chem
from rdkit.Chem import AllChem, Descriptors, QED

from .config import RANDOM_SEED

logger = logging.getLogger(__name__)


def smiles_to_mol3d(smiles: str) -> Chem.Mol:
    """SMILES → 3D RDKit molecule (explicit Hs, ETKDGv3 embed, MMFF minimised)."""
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        raise ValueError(f"Invalid SMILES: {smiles!r}")
    mol = Chem.AddHs(mol)

    params = AllChem.ETKDGv3()
    params.randomSeed = RANDOM_SEED
    if AllChem.EmbedMolecule(mol, params) == -1:
        params.useRandomCoords = True
        if AllChem.EmbedMolecule(mol, params) == -1:
            raise RuntimeError("Could not generate 3D coordinates for this molecule.")

    try:
        if AllChem.MMFFOptimizeMolecule(mol, maxIters=2000) == -1:
            AllChem.UFFOptimizeMolecule(mol, maxIters=2000)
    except Exception:
        pass
    return mol


def mol_to_pdbqt_string(mol: Chem.Mol) -> str:
    """3D RDKit mol (with Hs) → PDBQT string via Meeko, with finite partial charges."""
    try:
        from meeko import MoleculePreparation, PDBQTWriterLegacy
    except ImportError as exc:
        raise RuntimeError("Meeko is not installed (pip install meeko).") from exc

    # Large, highly charged and phosphorus-containing molecules can make
    # RDKit's Gasteiger solver return NaN/Inf. Supply checked charges through
    # Meeko's supported "read" model to prevent another implicit calculation.
    charges = None
    try:
        AllChem.ComputeGasteigerCharges(mol, nIter=12, throwOnParamFailure=False)
        values = [float(atom.GetProp("_GasteigerCharge")) for atom in mol.GetAtoms()]
        if all(math.isfinite(value) for value in values):
            charges = values
    except Exception:
        pass

    if charges is None:
        try:
            mmff = AllChem.MMFFGetMoleculeProperties(mol, mmffVariant="MMFF94")
            values = [float(mmff.GetMMFFPartialCharge(i)) for i in range(mol.GetNumAtoms())] if mmff else []
            if len(values) == mol.GetNumAtoms() and all(math.isfinite(value) for value in values):
                charges = values
        except Exception:
            pass

    if charges is None:
        logger.warning("Gasteiger and MMFF charges unavailable; using finite zero partial charges for PDBQT")
        charges = [0.0] * mol.GetNumAtoms()
    for atom, charge in zip(mol.GetAtoms(), charges):
        atom.SetDoubleProp("_TriposPartialCharge", charge)

    prep = MoleculePreparation(charge_model="read", charge_atom_prop="_TriposPartialCharge")
    setups = prep.prepare(mol)
    setup = setups[0] if isinstance(setups, (list, tuple)) else setups
    out = PDBQTWriterLegacy.write_string(setup)
    if isinstance(out, tuple):
        pdbqt, ok, err = out[0], out[1], (out[2] if len(out) > 2 else "")
        if not ok:
            raise RuntimeError(f"Meeko PDBQT write error: {err}")
        return pdbqt
    return out


def mol_to_sdf_block(mol: Chem.Mol) -> str:
    """3D RDKit mol → SDF/MOL block string (what Gnina & KarmaDock want)."""
    return Chem.MolToMolBlock(mol)


def compute_qed(mol: Chem.Mol) -> float:
    try:
        return round(QED.qed(Chem.RemoveHs(mol)), 4)
    except Exception:
        return 0.0


def mol_properties(mol: Chem.Mol) -> Dict[str, Any]:
    m = Chem.RemoveHs(mol)
    return {
        "mw": round(Descriptors.MolWt(m), 2),
        "logp": round(Descriptors.MolLogP(m), 2),
        "hbd": Descriptors.NumHDonors(m),
        "hba": Descriptors.NumHAcceptors(m),
        "tpsa": round(Descriptors.TPSA(m), 2),
        "rot_bonds": Descriptors.NumRotatableBonds(m),
    }
