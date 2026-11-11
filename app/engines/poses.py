# ══ CELL 12 · app/engines/poses.py ══
"""
poses.py — Convert docked poses (PDBQT / SDF / PDB) into clean PDB text.

The 3D viewer only needs PDB. Every engine returns its pose in a different
format, so everything is normalised here:

    pose_to_pdb(text)      -> PDB string (with bonds when possible)
    pose_centroid(text)    -> (x, y, z) of the ligand
    first_mol_from_sdf(t)  -> RDKit Mol (with Hs) or None
"""
import logging
import re
from typing import Optional, Tuple

from rdkit import Chem

logger = logging.getLogger(__name__)

# AutoDock atom types → chemical elements
_AD_TYPES = {
    "A": "C", "C": "C", "N": "N", "NA": "N", "NS": "N", "OA": "O", "OS": "O", "O": "O",
    "S": "S", "SA": "S", "HD": "H", "HS": "H", "H": "H", "P": "P", "F": "F",
    "CL": "Cl", "BR": "Br", "I": "I", "ZN": "Zn", "FE": "Fe", "MG": "Mg", "CA": "Ca",
    "MN": "Mn", "CU": "Cu", "NI": "Ni", "CO": "Co", "CD": "Cd", "B": "B", "SE": "Se",
}


def detect_format(text: str) -> str:
    """Return 'sdf', 'pdbqt' or 'pdb'."""
    if "$$$$" in text or "V2000" in text or "V3000" in text:
        return "sdf"
    for line in text.splitlines():
        if line.startswith(("ROOT", "ENDROOT", "BRANCH", "TORSDOF")):
            return "pdbqt"
    for line in text.splitlines():
        if line.startswith(("ATOM", "HETATM")) and len(line) >= 78:
            tail = line[66:].split()
            if len(tail) >= 2 and re.match(r"^[+-]?\d", tail[0]):
                return "pdbqt"
    return "pdb"


def first_mol_from_sdf(text: str) -> Optional[Chem.Mol]:
    """First molecule of an SDF/MOL string (keeps Hs). Tolerates a missing title line."""
    blocks = text.split("$$$$")
    raw = blocks[0] if blocks else text
    variants = [raw, raw.lstrip("\n"), "\n" + raw.lstrip("\n")]
    for block in variants:
        for sanitize in (True, False):
            try:
                mol = Chem.MolFromMolBlock(block, sanitize=sanitize, removeHs=False)
                if mol is not None and mol.GetNumConformers():
                    return mol
            except Exception:
                continue
    return None


def _first_model(text: str) -> str:
    """Return only the first MODEL…ENDMDL block (or the whole text)."""
    out, in_model, seen = [], False, False
    for line in text.splitlines():
        if line.startswith("MODEL"):
            if seen:
                break
            in_model = seen = True
            continue
        if line.startswith("ENDMDL"):
            break
        out.append(line)
    return "\n".join(out)


def _pdbqt_to_pdb_manual(text: str) -> str:
    """Fallback PDBQT→PDB: fixed columns, element from AutoDock type (no bonds)."""
    lines, serial = [], 0
    for line in _first_model(text).splitlines():
        if not line.startswith(("ATOM", "HETATM")):
            continue
        try:
            x, y, z = float(line[30:38]), float(line[38:46]), float(line[46:54])
        except ValueError:
            continue
        serial += 1
        name = line[12:16]
        resn = line[17:20].strip() or "LIG"
        chain = line[21] if len(line) > 21 and line[21].strip() else "A"
        try:
            resi = int(line[22:26])
        except ValueError:
            resi = 1
        atype = line[77:79].strip().upper() if len(line) >= 78 else ""
        elem = _AD_TYPES.get(atype) or re.sub(r"[^A-Za-z]", "", name)[:1].upper() or "C"
        lines.append(
            f"HETATM{serial:5d} {name:<4s} {resn:>3s} {chain}{resi:4d}    "
            f"{x:8.3f}{y:8.3f}{z:8.3f}  1.00  0.00          {elem:>2s}"
        )
    lines.append("END")
    return "\n".join(lines) + "\n"


def pdbqt_to_mol(text: str) -> Optional[Chem.Mol]:
    """PDBQT (as written by Meeko/Vina) → RDKit mol with bonds, via Meeko."""
    try:
        from meeko import PDBQTMolecule, RDKitMolCreate
        pm = PDBQTMolecule(_first_model_wrapped(text), skip_typing=True)
        mols = RDKitMolCreate.from_pdbqt_mol(pm)
        mols = [m for m in mols if m is not None]
        return mols[0] if mols else None
    except Exception as exc:
        logger.debug("Meeko PDBQT→mol failed: %s", exc)
        return None


def _first_model_wrapped(text: str) -> str:
    """Keep MODEL 1 with its REMARK lines so Meeko can read it."""
    if "MODEL" not in text:
        return text
    out, seen_model = [], 0
    for line in text.splitlines():
        if line.startswith("MODEL"):
            seen_model += 1
            if seen_model > 1:
                break
        out.append(line)
        if line.startswith("ENDMDL"):
            break
    return "\n".join(out) + "\n"


def pose_to_pdb(text: str) -> str:
    """Any pose text → PDB text suitable for 3Dmol. Never raises."""
    if not text or not text.strip():
        return ""
    try:
        fmt = detect_format(text)
        if fmt == "sdf":
            mol = first_mol_from_sdf(text)
            if mol is not None:
                return Chem.MolToPDBBlock(mol)
        elif fmt == "pdbqt":
            mol = pdbqt_to_mol(text)
            if mol is not None:
                return Chem.MolToPDBBlock(mol)
            return _pdbqt_to_pdb_manual(text)
        else:
            return _first_model(text) + "\nEND\n"
    except Exception as exc:
        logger.warning("pose_to_pdb failed (%s) — trying manual conversion", exc)
        try:
            return _pdbqt_to_pdb_manual(text)
        except Exception:
            pass
    return ""


def pose_coords(text: str):
    """Heavy-atom coordinates of the first pose as a list of (x, y, z)."""
    pdb = pose_to_pdb(text)
    pts = []
    for line in pdb.splitlines():
        if line.startswith(("ATOM", "HETATM")):
            elem = line[76:78].strip().upper()
            if elem == "H":
                continue
            try:
                pts.append((float(line[30:38]), float(line[38:46]), float(line[46:54])))
            except ValueError:
                pass
    return pts


def pose_centroid(text: str) -> Optional[Tuple[float, float, float]]:
    pts = pose_coords(text)
    if not pts:
        return None
    n = len(pts)
    return (sum(p[0] for p in pts) / n, sum(p[1] for p in pts) / n, sum(p[2] for p in pts) / n)
