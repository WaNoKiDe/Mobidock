"""
poses.py — compare docking poses (RMSD) and find the experimental ligand pose of a PDB entry.

All poses live in the receptor's coordinate frame, so RMSD is measured WITHOUT re-aligning them.
Molecules are read as heavy atoms with bonds guessed from distances; the common core of two molecules is found
with RDKit's maximum-common-substructure search, so a modified molecule can be compared with its parent.
"""
import logging
import math
import re
import time
from typing import Optional, Tuple

from rdkit import Chem
from rdkit.Chem import rdFMCS

from . import recorddb

logger = logging.getLogger("medchem.poses")
_NOT_LIGANDS = set("""HOH GOL EDO PEG PGE PG4 1PE P6G SO4 PO4 ACT ACY FMT DMS TRS MPD BME IPA EOH MES CIT TLA FLC NO3 CL NA K CA MG ZN
MN FE CU CO NI CD HG IOD BR NAG NDG MAN BMA FUC GAL GLC SIA BGC SGN CO3 OXY PEO NH4 SCN AZI HEM FAD FMN NAD NAP NDP NAI COA SAM
SAH PLP TPP BOG OLC DTT EPE PGO PE4 PE5 PE8 BU3 IMD MLI MLA TAR SUC TRE CAC""".split())
_FAIL: dict = {}


def pose_mol(pdb_text: str):
    """Heavy-atom molecule with the coordinates of a pose (PDB text), or None."""
    lines = [l for l in (pdb_text or "").splitlines() if l.startswith(("ATOM", "HETATM", "CONECT"))]
    if not any(l.startswith(("ATOM", "HETATM")) for l in lines):
        return None
    try:
        m = Chem.MolFromPDBBlock("\n".join(lines) + "\nEND\n", sanitize=False, removeHs=False, proximityBonding=True)
        if m is None or m.GetNumConformers() == 0:
            return None
        m = Chem.RemoveHs(m, sanitize=False)
        m.UpdatePropertyCache(strict=False)
        Chem.FastFindRings(m)
        return m if m.GetNumAtoms() else None
    except Exception:
        return None


def rmsd_common(a, b, timeout: int = 5, min_atoms: int = 6) -> Optional[Tuple[float, int]]:
    """(RMSD in A, atoms compared) over the largest common core of two poses, without re-aligning; None if too little in common."""
    if a is None or b is None:
        return None
    try:
        r = rdFMCS.FindMCS([a, b], timeout=timeout, atomCompare=rdFMCS.AtomCompare.CompareElements,
                           bondCompare=rdFMCS.BondCompare.CompareAny, ringMatchesRingOnly=True, completeRingsOnly=True)
    except Exception:
        return None
    if r.numAtoms < min_atoms or r.numAtoms < 0.5 * min(a.GetNumAtoms(), b.GetNumAtoms()):
        return None
    q = Chem.MolFromSmarts(r.smartsString)
    if q is None:
        return None
    ma = a.GetSubstructMatches(q, uniquify=False, maxMatches=200)
    mb = b.GetSubstructMatches(q, uniquify=False, maxMatches=200)
    if not ma or not mb:
        return None
    pa, pb = a.GetConformer().GetPositions(), b.GetConformer().GetPositions()
    best = None
    for x in ma[:40]:                                  # molecules with symmetry match several ways: keep the best correspondence
        for y in mb[:40]:
            d = pa[list(x)] - pb[list(y)]
            v = math.sqrt(float((d * d).sum()) / len(x))
            best = v if best is None or v < best else best
    return (best, int(r.numAtoms)) if best is not None else None


def pdb_for_target(target_id: str) -> Optional[str]:
    try:
        from .receptor_prep import RECEPTORS
        cand = str((RECEPTORS.get(target_id) or {}).get("pdb_id", ""))
    except Exception:
        cand = ""
    for c in (cand, target_id):
        if re.fullmatch(r"[0-9][A-Za-z0-9]{3}", c or ""):
            return c.upper()
    m = re.match(r"PDB_([0-9][A-Za-z0-9]{3})_", target_id or "")
    return m.group(1).upper() if m else None


def native_ligand_pdb(raw_pdb: str, near=None) -> Tuple[str, str]:
    """(PDB lines, residue name) of the drug-like ligand in a PDB file: the one closest to `near` (x, y, z) or the largest."""
    groups: dict = {}
    for ln in raw_pdb.splitlines():
        if not ln.startswith("HETATM"):
            continue
        res = ln[17:20].strip()
        if res in _NOT_LIGANDS:
            continue
        el = (ln[76:78].strip() or ln[12:16].strip()[:1]).upper()
        if el == "H":
            continue
        groups.setdefault((res, ln[21], ln[22:27]), []).append(ln)
    best, best_key = None, None
    for key, ls in groups.items():
        if len(ls) < 8:
            continue
        try:
            xs = [float(l[30:38]) for l in ls]; ys = [float(l[38:46]) for l in ls]; zs = [float(l[46:54]) for l in ls]
        except ValueError:
            continue
        cen = (sum(xs) / len(xs), sum(ys) / len(ys), sum(zs) / len(zs))
        score = math.dist(cen, near) if near else -len(ls)
        if best is None or score < best:
            best, best_key = score, key
    if best_key is None:
        return "", ""
    return "\n".join(groups[best_key]), best_key[0]


def native_mol(target_id: str):
    """The experimental ligand pose of the target's PDB entry (downloaded once, then kept in the database), or None."""
    row = recorddb.get_native(target_id)
    if row:
        return pose_mol(row["pdb"])
    if time.time() - _FAIL.get(target_id, 0) < 600:
        return None
    pdb = pdb_for_target(target_id)
    if not pdb:
        return None
    try:
        from . import receptor_prep as rp
        raw = rp.download_pdb(pdb)
        near = (getattr(rp, "KNOWN_SITES", {}).get(target_id) or {}).get("center")
        lines, comp = native_ligand_pdb(raw or "", near)
        if not lines:
            _FAIL[target_id] = time.time()
            return None
        recorddb.save_native(target_id, comp, lines)
        return pose_mol(lines)
    except Exception as exc:
        logger.info("experimental ligand pose for %s unavailable: %s", target_id, exc)
        _FAIL[target_id] = time.time()
        return None
