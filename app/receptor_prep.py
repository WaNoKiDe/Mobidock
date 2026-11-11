# ══ CELL 08 · app/receptor_prep.py ══
"""
receptor_prep.py — Download & prepare protein receptors for docking.

For every target this writes three files into receptors/:
    <name>_clean.pdb   heavy atoms only (shown in the 3D viewer)
    <name>_H.pdb       protonated at pH 7.4 (used by Gnina / DiffDock)
    <name>.pdbqt       AutoDock-typed + Gasteiger charges (used by Vina)

Conversion uses Open Babel (pip install openbabel-wheel). There is no
"copy the PDB and hope" fallback any more — if conversion fails you get an error.

Imported / uploaded proteins also get an automatic search box centred on the
largest co-crystallised ligand (falls back to the protein centre).
"""
import logging
import subprocess
import tempfile
import threading
from pathlib import Path
from typing import Dict, List, Optional, Tuple

from .config import RECEPTORS_DIR, DEFAULT_BOX_SIZE

logger = logging.getLogger(__name__)

# Targets are created only when the user imports or uploads a PDB structure.
RECEPTORS = {}

# Residue names that are NOT drug-like ligands (water, ions, buffers, cryoprotectants…)
_NOT_LIGANDS = {
    "HOH", "WAT", "DOD", "NA", "CL", "K", "MG", "ZN", "CA", "MN", "FE", "FE2", "CU", "NI", "CO", "CD",
    "SO4", "PO4", "GOL", "EDO", "PEG", "PGE", "DMS", "ACT", "FMT", "TRS", "MPD", "BME", "IOD", "BR",
    "NO3", "NH4", "ACE", "NH2", "PG4", "1PE", "P6G", "EPE", "MES", "CIT", "TLA", "SCN", "AZI", "IMD",
}


# ── Download & clean ──────────────────────────────────────────────────────────

def download_pdb(pdb_id: str) -> Optional[str]:
    import urllib.request
    url = f"https://files.rcsb.org/download/{pdb_id.upper()}.pdb"
    try:
        req = urllib.request.Request(url, headers={"User-Agent": "MedChemGame/3.1"})
        with urllib.request.urlopen(req, timeout=40) as resp:
            return resp.read().decode("utf-8", errors="ignore")
    except Exception as exc:
        logger.error("Failed to download %s: %s", pdb_id, exc)
        return None


def clean_pdb(pdb_text: str) -> str:
    """Keep protein ATOM records (first model, altloc A/blank); drop water, ligands, ions."""
    lines = []
    for line in pdb_text.splitlines():
        if line.startswith("ENDMDL"):
            break                              # only the first model (NMR structures)
        if line.startswith("ATOM"):
            alt = line[16] if len(line) > 16 else " "
            if alt in (" ", "A", "1"):
                lines.append(line[:16] + " " + line[17:])
        elif line.startswith("TER"):
            lines.append("TER")
    lines.append("END")
    return "\n".join(lines) + "\n"


# ── Binding-site detection ────────────────────────────────────────────────────

def find_binding_site(raw_pdb: str) -> Dict:
    """
    Pick the search box automatically.
    1. Largest co-crystallised drug-like ligand → centre on it, box = ligand extent + 10 Å.
    2. Otherwise → centre of the protein, 30 Å box (user should click to refine).
    """
    groups: Dict[tuple, List[tuple]] = {}
    protein_pts: List[tuple] = []
    for line in raw_pdb.splitlines():
        rec = line[:6]
        if rec not in ("ATOM  ", "HETATM"):
            continue
        try:
            x, y, z = float(line[30:38]), float(line[38:46]), float(line[46:54])
        except ValueError:
            continue
        if rec == "ATOM  ":
            protein_pts.append((x, y, z))
            continue
        resn = line[17:20].strip().upper()
        elem = line[76:78].strip().upper() if len(line) >= 78 else ""
        if resn in _NOT_LIGANDS or elem == "H":
            continue
        groups.setdefault((resn, line[21], line[22:26].strip()), []).append((x, y, z))

    ligs = {k: v for k, v in groups.items() if len(v) >= 8}
    if ligs:
        key = max(ligs, key=lambda k: len(ligs[k]))
        pts = ligs[key]
        n = len(pts)
        c = tuple(round(sum(p[i] for p in pts) / n, 2) for i in range(3))
        ext = [max(p[i] for p in pts) - min(p[i] for p in pts) for i in range(3)]
        size = tuple(round(min(40.0, max(18.0, e + 10.0)), 1) for e in ext)
        return {"center": c, "size": size, "source": "co-crystal ligand", "ligand": key[0]}

    if protein_pts:
        n = len(protein_pts)
        c = tuple(round(sum(p[i] for p in protein_pts) / n, 2) for i in range(3))
        return {"center": c, "size": (30.0, 30.0, 30.0), "source": "protein centre", "ligand": None}
    return {"center": (0.0, 0.0, 0.0), "size": DEFAULT_BOX_SIZE, "source": "default", "ligand": None}


# ── PDB → protonated PDB + PDBQT (Open Babel) ─────────────────────────────────

def is_valid_pdbqt(path: Path) -> bool:
    """A real receptor PDBQT has AutoDock atom types in columns 78-79."""
    try:
        n_ok = n_atoms = 0
        with open(path, "r", encoding="utf-8", errors="ignore") as fh:
            for line in fh:
                if line.startswith(("ATOM", "HETATM")):
                    n_atoms += 1
                    if len(line) >= 78 and line[77:79].strip() and line[70:76].strip():
                        n_ok += 1
                    if n_atoms >= 50:
                        break
        return n_atoms > 0 and n_ok == n_atoms
    except Exception:
        return False


_OB_LOCK = threading.Lock()      # Open Babel's SWIG bindings are not thread-safe


def _prep_with_pybel(clean_path: Path, h_pdb: Path, pdbqt: Path) -> None:
    """Protonate at pH 7.4 and write PDB + rigid-receptor PDBQT using the OBConversion API."""
    from openbabel import openbabel as ob
    with _OB_LOCK:
        ob.obErrorLog.SetOutputLevel(0)
        mol = ob.OBMol()
        conv = ob.OBConversion()
        conv.SetInFormat("pdb")
        if not conv.ReadFile(mol, str(clean_path)):
            raise RuntimeError("Open Babel could not read the PDB file")
        mol.AddHydrogens(False, True, 7.4)           # polaronly=False, correctForPH=True, pH 7.4

        conv.SetOutFormat("pdb")
        conv.WriteFile(mol, str(h_pdb))

        out = ob.OBConversion()
        out.SetOutFormat("pdbqt")
        out.AddOption("r", ob.OBConversion.OUTOPTIONS)   # rigid receptor (no ROOT/BRANCH)
        out.WriteFile(mol, str(pdbqt))


def _prep_with_obabel_cli(clean_path: Path, h_pdb: Path, pdbqt: Path) -> None:
    subprocess.run(["obabel", str(clean_path), "-O", str(h_pdb), "-h", "-p", "7.4"],
                   check=True, capture_output=True, timeout=300)
    subprocess.run(["obabel", str(h_pdb), "-O", str(pdbqt), "-xr"],
                   check=True, capture_output=True, timeout=300)


def convert_receptor(clean_path: Path, h_pdb: Path, pdbqt: Path) -> None:
    """Write h_pdb and pdbqt from clean_path. Raises RuntimeError on failure."""
    errors = []
    for fn in (_prep_with_pybel, _prep_with_obabel_cli):
        try:
            fn(clean_path, h_pdb, pdbqt)
            if is_valid_pdbqt(pdbqt):
                return
            errors.append(f"{fn.__name__}: output was not valid PDBQT")
        except Exception as exc:
            errors.append(f"{fn.__name__}: {exc}")
    raise RuntimeError(
        "Receptor conversion failed. Install Open Babel (pip install openbabel-wheel). Details: "
        + " | ".join(errors)
    )


def _paths(stem: str):
    return (RECEPTORS_DIR / f"{stem}_clean.pdb", RECEPTORS_DIR / f"{stem}_H.pdb", RECEPTORS_DIR / f"{stem}.pdbqt")


# ── Target preparation ────────────────────────────────────────────────────────

def prepare_receptor(target_id: str, force: bool = False) -> Optional[Path]:
    """Prepare a registered target. Returns the PDBQT path or None."""
    RECEPTORS_DIR.mkdir(parents=True, exist_ok=True)
    info = RECEPTORS[target_id]
    stem = f"{target_id}_{info['pdb_id']}"
    clean, h_pdb, pdbqt = _paths(stem)

    if not force and pdbqt.exists() and h_pdb.exists() and clean.exists() and is_valid_pdbqt(pdbqt):
        return pdbqt
    if pdbqt.exists() and not is_valid_pdbqt(pdbqt):
        logger.warning("%s: old PDBQT is invalid (plain PDB) — rebuilding", pdbqt.name)

    raw = download_pdb(info["pdb_id"])
    if not raw:
        logger.error("Could not download %s", info["pdb_id"])
        return None
    clean.write_text(clean_pdb(raw), encoding="utf-8")
    try:
        convert_receptor(clean, h_pdb, pdbqt)
    except Exception as exc:
        logger.error("%s: %s", target_id, exc)
        return None
    logger.info("Receptor %s ready: %s", target_id, pdbqt.name)
    return pdbqt


def prepare_all_receptors() -> dict:
    return {tid: prepare_receptor(tid) for tid in RECEPTORS}


# ── Custom / imported targets ─────────────────────────────────────────────────

def process_custom_pdb(target_id: str, raw_pdb: str, ligand_comp_id: Optional[str] = None) -> dict:
    """Prepare an uploaded or RCSB-imported PDB; auto-pick the search box."""
    RECEPTORS_DIR.mkdir(parents=True, exist_ok=True)
    clean, h_pdb, pdbqt = _paths(target_id)
    raw_path = RECEPTORS_DIR / f"{target_id}_source.pdb"
    raw_path.write_text(raw_pdb, encoding="utf-8")
    cleaned = clean_pdb(raw_pdb)
    if cleaned.count("ATOM") == 0:
        raise ValueError("No protein ATOM records found in this PDB file.")
    clean.write_text(cleaned, encoding="utf-8")

    site = (
        find_component_binding_site(raw_pdb, ligand_comp_id)
        if ligand_comp_id else find_binding_site(raw_pdb)
    )
    convert_receptor(clean, h_pdb, pdbqt)

    return {
        "pdb_id": target_id,
        "name": "Custom Target",
        "disease": "User Upload",
        "center": site["center"],
        "size": site["size"],
        "site_source": site["source"],
        "site_ligand": site["ligand"],
        "color": "#1abc9c",
        "icon": "🧪",
        "pdbqt_path": pdbqt,
        "pdb_path": clean,
        "raw_pdb_path": raw_path,
    }


# ── CLI ───────────────────────────────────────────────────────────────────────
if __name__ == "__main__":
    logging.basicConfig(level=logging.INFO, format="%(levelname)s: %(message)s")
    for tid, path in prepare_all_receptors().items():
        print(f"  [{'OK' if path and path.exists() else 'FAIL'}] {tid}: {path}")

# ==== 49t known ligand sites begin ====
# The site of the PDB entry's own ligand, remembered when the target is loaded. Moving the grid box later changes
# RECEPTORS[...]["center"], but never this.
KNOWN_SITES = {k: {"center": tuple(v["center"]), "size": tuple(v["size"])} for k, v in RECEPTORS.items()}
_orig_process_custom_pdb = process_custom_pdb


def process_custom_pdb(*args, **kwargs):
    info = _orig_process_custom_pdb(*args, **kwargs)
    try:
        tid = args[0] if args else kwargs.get("target_id")
        if tid and info.get("site_source") == "co-crystal ligand":
            KNOWN_SITES[tid] = {"center": tuple(info["center"]), "size": tuple(info["size"])}
    except Exception:
        pass
    return info
# ==== 49t known ligand sites end ====


# === LIGAND_PDB_COMPONENT_SITE_V1 ===
def find_component_binding_site(raw_pdb: str, ligand_comp_id: str, chain=None, residue=None, insertion=None) -> Dict:
    """
    Centre the docking grid on the actual co-crystallized ligand selected from RCSB.
    The receptor is still cleaned separately, so this ligand is removed before docking.
    """
    comp_id = (ligand_comp_id or "").strip().upper()
    if not comp_id:
        raise ValueError("The selected PDB result has no chemical-component ID.")

    groups: Dict[tuple, List[tuple]] = {}
    for line in raw_pdb.splitlines():
        if not line.startswith("HETATM"):
            continue

        resn = line[17:20].strip().upper()
        if resn != comp_id:
            continue

        element = line[76:78].strip().upper() if len(line) >= 78 else ""
        if element == "H":
            continue

        try:
            xyz = (float(line[30:38]), float(line[38:46]), float(line[46:54]))
        except ValueError:
            continue

        atom_chain = line[21].strip() or "_"
        atom_residue = line[22:26].strip()
        atom_insertion = line[26].strip()
        if chain is not None and atom_chain != str(chain):
            continue
        if residue is not None and atom_residue != str(residue):
            continue
        if insertion is not None and atom_insertion != str(insertion or ""):
            continue
        key = (atom_chain, atom_residue, atom_insertion, resn)
        groups.setdefault(key, []).append(xyz)

    if not groups:
        raise ValueError(
            f"Could not locate co-crystallized ligand {comp_id!r} in the PDB coordinate file."
        )

    # If a crystal has several copies, use the best-defined / largest copy.
    selected_key, points = max(groups.items(), key=lambda item: len(item[1]))
    n = len(points)
    center = tuple(round(sum(p[i] for p in points) / n, 2) for i in range(3))
    extent = [max(p[i] for p in points) - min(p[i] for p in points) for i in range(3)]

    # Native ligand extent plus 12 Å padding, within Vina's laptop-safe box range.
    size = tuple(round(min(30.0, max(18.0, edge + 12.0)), 1) for edge in extent)

    return {
        "center": center,
        "size": size,
        "source": f"selected co-crystal ligand {comp_id} at chain {selected_key[0]} residue {selected_key[1]}",
        "ligand": comp_id,
        "chain": selected_key[0],
        "residue": selected_key[1],
        "insertion": selected_key[2],
    }
# === END LIGAND_PDB_COMPONENT_SITE_V1 ===
