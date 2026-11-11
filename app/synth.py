"""
synth.py — Synthesizability: "could a chemist actually make this molecule?"

Method: the SA score of Ertl & Schuffenhauer (J. Cheminformatics 2009), shipped with RDKit (Contrib/SA_Score).
    1  = trivially easy to make          10 = practically impossible
It combines (a) how common the molecule's ring/substituent patterns are among known compounds with
(b) penalties for size, stereocentres, spiro/bridged rings and macrocycles.

This is a fast, 2-D HEURISTIC, not a retrosynthesis. A low score does not guarantee a route exists, and
some real drugs (natural products especially) score high. It is meant to steer players away from
"great docking score, impossible molecule" answers.

Modes for the gate (same idea as the Lipinski gate):
    standard = reject only "very hard" molecules (SA > 6.5)
    strict   = reject anything harder than "moderate" (SA > 5.0)
    off      = report only, never reject
"""
import importlib.util
import logging
import math
import sys
import threading
import urllib.request
from pathlib import Path
from typing import Any, Dict, Optional

from rdkit import Chem
from rdkit.Chem import rdMolDescriptors

from .config import DATA_DIR

logger = logging.getLogger(__name__)

SYNTH_MODES = ("standard", "strict", "off")
SA_LIMITS = {"standard": 6.5, "strict": 5.0}
SYNTH_WEIGHT = 2.0          # game-score points for a perfectly easy molecule (QED is worth up to 3.0)

# (upper SA bound, key, label, colour)
TIERS = (
    (3.0, "easy",      "🟢 Easy to make",   "#2af598"),
    (5.0, "moderate",  "🟡 Moderate",       "#ffd60a"),
    (6.5, "hard",      "🟠 Challenging",    "#f8961e"),
    (99., "very_hard", "🔴 Very hard",      "#ff4d6d"),
)
_K = 9.0 / 6.5              # how many SA points one unit of "raw" penalty is worth (from the SA formula)

_SA_URLS = {
    "sascorer.py": "https://raw.githubusercontent.com/rdkit/rdkit/master/Contrib/SA_Score/sascorer.py",
    "fpscores.pkl.gz": "https://raw.githubusercontent.com/rdkit/rdkit/master/Contrib/SA_Score/fpscores.pkl.gz",
}
_sa_module = None
_sa_error: Optional[str] = None
_lock = threading.Lock()


# ── loading the SA scorer (RDKit Contrib → RDConfig path → one-time download) ──
def _load_from_file(path: Path):
    spec = importlib.util.spec_from_file_location("medchem_sascorer", str(path))
    if spec is None or spec.loader is None:
        raise ImportError(f"Could not load SA scorer from {path}")
    mod = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(mod)
    return mod


def _prime_sascorer(module, score_file: Optional[Path] = None):
    """Load the fragment table now, so a missing data file cannot fail on first scoring."""
    if score_file is not None:
        module.readFragmentScores(str(score_file))
    elif getattr(module, "_fscores", None) is None:
        module.readFragmentScores()
    scores = getattr(module, "_fscores", None)
    if not isinstance(scores, dict) or not scores:
        raise RuntimeError("RDKit SA fragment score table is empty")
    return module


def _download_sa_assets(directory: Path, force: bool = False):
    directory.mkdir(parents=True, exist_ok=True)
    for name, url in _SA_URLS.items():
        path = directory / name
        if force or not path.exists():
            logger.info("Downloading %s …", name)
            with urllib.request.urlopen(url, timeout=60) as response:
                path.write_bytes(response.read())


def _load_sascorer():
    global _sa_module, _sa_error
    if _sa_module is not None:
        return _sa_module
    with _lock:
        if _sa_module is not None:
            return _sa_module
        failures = []
        try:
            from rdkit.Contrib.SA_Score import sascorer
            _sa_module = _prime_sascorer(sascorer)
            _sa_error = None
            return _sa_module
        except Exception as exc:
            failures.append(exc)

        try:
            from rdkit import RDConfig
            scorer_path = Path(RDConfig.RDContribDir) / "SA_Score" / "sascorer.py"
            if scorer_path.exists():
                _sa_module = _prime_sascorer(_load_from_file(scorer_path))
                _sa_error = None
                return _sa_module
        except Exception as exc:
            failures.append(exc)

        # Some RDKit installs ship sascorer.py without fpscores.pkl.gz. Keep a
        # writable app-local copy so the score can still be calculated there.
        local_dir = DATA_DIR / "sa_score"
        try:
            _download_sa_assets(local_dir)
            local_score_file = local_dir / "fpscores.pkl.gz"
            try:
                _sa_module = _prime_sascorer(_load_from_file(local_dir / "sascorer.py"), local_score_file)
            except Exception:
                # Repair partial or stale downloads once before reporting failure.
                _download_sa_assets(local_dir, force=True)
                _sa_module = _prime_sascorer(_load_from_file(local_dir / "sascorer.py"), local_score_file)
            _sa_error = None
        except Exception as exc:
            failures.append(exc)
            detail = failures[-1]
            _sa_error = f"SA scorer unavailable ({type(detail).__name__}: {detail})"
            logger.warning("%s", _sa_error)
        return _sa_module


def sa_available() -> bool:
    return _load_sascorer() is not None


# ── helpers ──
def _tier(sa: float):
    for upper, key, label, colour in TIERS:
        if sa <= upper:
            return key, label, colour
    return TIERS[-1][1:]


def _to_mol(x) -> Optional[Chem.Mol]:
    if isinstance(x, Chem.Mol):
        m = x
    else:
        s = (x or "").strip()
        m = Chem.MolFromSmiles(s) if s else None
    if m is None:
        return None
    return Chem.RemoveHs(m)


def _breakdown(mol: Chem.Mol, sa: float) -> list:
    """Which structural features push the SA score up, in SA points (same formula the score uses)."""
    n_atoms = mol.GetNumAtoms()
    n_stereo = len(Chem.FindMolChiralCenters(mol, includeUnassigned=True))
    n_spiro = rdMolDescriptors.CalcNumSpiroAtoms(mol)
    n_bridge = rdMolDescriptors.CalcNumBridgeheadAtoms(mol)
    n_macro = sum(1 for r in mol.GetRingInfo().AtomRings() if len(r) > 8)
    parts = {
        "size":   _K * (n_atoms ** 1.005 - n_atoms),
        "stereo": _K * math.log10(n_stereo + 1),
        "spiro":  _K * math.log10(n_spiro + 1),
        "bridge": _K * math.log10(n_bridge + 1),
        "macro":  _K * math.log10(2) if n_macro else 0.0,
    }
    text = {
        "size":   ("📏", "Size", f"{n_atoms} heavy atoms: every extra atom is more bonds to make. Trim what the pocket does not need."),
        "stereo": ("🪞", "Stereocentres", f"{n_stereo} stereocentre(s): each one must be controlled in the lab. Remove or simplify some."),
        "spiro":  ("🌀", "Spiro rings", f"{n_spiro} spiro atom(s): rings sharing one atom are awkward to build."),
        "bridge": ("🌉", "Bridged rings", f"{n_bridge} bridgehead atom(s): cage-like ring systems are hard to assemble."),
        "macro":  ("⭕", "Macrocycle", "A ring of more than 8 atoms: ring-closing big rings is notoriously low-yielding."),
    }
    out = []
    for key in ("size", "stereo", "spiro", "bridge", "macro"):
        pts = parts[key]
        if pts >= 0.15:
            icon, label, tip = text[key]
            out.append({"key": key, "icon": icon, "label": label, "points": round(pts, 2), "tip": tip})
    rest = sa - sum(parts.values())
    out.append({"key": "fragments", "icon": "🧩", "label": "Building blocks",
                "points": round(rest, 2),
                "tip": "How common the ring and substituent patterns are among known compounds (rare patterns = harder), "
                       "after a small bonus for symmetry."})
    out.sort(key=lambda d: d["points"], reverse=True)
    return out


# Explanation uses the same sparse Morgan IDs and fragment table as RDKit SA.
def _fragment_explanation(mol, sas):
    from rdkit.Chem import rdFingerprintGenerator, AllChem, Draw
    import base64
    if sas._fscores is None:
        sas.readFragmentScores()
    if hasattr(sas, 'mfpgen'):
        extra = rdFingerprintGenerator.AdditionalOutput()
        extra.AllocateBitInfoMap()
        fp = sas.mfpgen.GetSparseCountFingerprint(mol, additionalOutput=extra)
        info = extra.GetBitInfoMap()
    else:
        info = {}
        fp = AllChem.GetMorganFingerprint(mol, 2, bitInfo=info)
    counts = fp.GetNonzeroElements()
    total = sum(counts.values())
    rows = []
    for bit, count in counts.items():
        weight = float(sas._fscores.get(bit, -4.0))
        rows.append({'id': str(bit), 'count': int(count),
                     'matched': bit in sas._fscores, 'weight': weight,
                     'contribution': count * weight / total,
                     'environments': [{'atom': int(a), 'radius': int(r)} for a,r in info.get(bit, ())]})
    rows.sort(key=lambda r: (-abs(r['contribution']), r['id']))
    # Show representative environments for the largest absolute contributions.
    # All IDs and all occurrences remain in the table, including unmapped IDs.
    for row in rows[:6]:
        try:
            svg = Draw.DrawMorganBit(mol, int(row['id']), info, useSVG=True)
            row['image'] = 'data:image/svg+xml;base64,' + base64.b64encode(svg.encode()).decode()
        except Exception:
            row['image'] = None
    n = mol.GetNumAtoms()
    penalties = {
        'Size': n ** 1.005 - n,
        'Stereocentres': math.log10(len(Chem.FindMolChiralCenters(mol, includeUnassigned=True)) + 1),
        'Spiro atoms': math.log10(rdMolDescriptors.CalcNumSpiroAtoms(mol) + 1),
        'Bridgehead atoms': math.log10(rdMolDescriptors.CalcNumBridgeheadAtoms(mol) + 1),
        'Macrocycles': math.log10(2) if any(len(r)>8 for r in mol.GetRingInfo().AtomRings()) else 0.0,
    }
    fragment = sum(r['contribution'] for r in rows)
    penalty = sum(penalties.values())
    density = 0.5 * math.log(n / len(counts)) if n > len(counts) else 0.0
    raw = fragment - penalty + density
    mapped = 11.0 - (raw + 5.0) * 9.0 / 6.5
    smoothed = 8.0 + math.log(mapped - 8.0) if mapped > 8.0 else mapped
    reconstructed = max(1.0, min(10.0, smoothed))
    actual = float(sas.calculateScore(mol))
    return {'rows': rows, 'unique': len(rows), 'occurrences': total,
            'matched_unique': sum(r['matched'] for r in rows),
            'matched_occurrences': sum(r['count'] for r in rows if r['matched']),
            'fragment_mean': fragment, 'penalties': penalties, 'penalty_total': penalty,
            'density_correction': density, 'raw': raw, 'mapped': mapped,
            'reconstructed_sa': reconstructed, 'verified': abs(reconstructed-actual)<1e-8}


def _safe_fragment_explanation(mol, sas):
    try:
        return _fragment_explanation(mol, sas)
    except Exception as exc:
        logger.warning('SA explanation unavailable: %s', exc)
        return {'error': 'Fragment explanation unavailable for this scorer version.'}

# ── public API ──
def analyze(smiles_or_mol, mode: str = "standard") -> Dict[str, Any]:
    """Full synthesizability report + gate verdict. Never raises for a missing scorer: check ['available']."""
    mode = (mode or "standard").strip().lower()
    if mode not in SYNTH_MODES:
        raise ValueError("Synthesizability mode must be 'standard', 'strict' or 'off'")
    mol = _to_mol(smiles_or_mol)
    if mol is None or mol.GetNumAtoms() == 0:
        raise ValueError("Cannot score an empty or invalid molecule")
    sa_smiles = Chem.MolToSmiles(mol, canonical=True)

    sas = _load_sascorer()
    if sas is None:
        return {"available": False, "error": _sa_error or "SA scorer unavailable", "mode": mode,
                "passed": True, "enforced": False, "skipped": True, "sa_smiles": sa_smiles}

    sa = round(float(sas.calculateScore(mol)), 2)
    norm = round((10.0 - sa) / 9.0, 4)                      # 1.0 = trivial, 0.0 = impossible
    key, label, colour = _tier(sa)
    limit = SA_LIMITS.get(mode)
    passed = True if limit is None else sa <= limit
    return {
        "available": True,
        "sa_smiles": sa_smiles,
        "method": "SA score (Ertl & Schuffenhauer 2009)",
        "sa": sa,
        "score": int(round(norm * 100)),                    # 0-100, higher = easier
        "norm": norm,
        "tier": key, "tier_label": label, "colour": colour,
        "drivers": _breakdown(mol, sa),
        "explanation": _safe_fragment_explanation(mol, sas),
        "mode": mode, "limit": limit, "enforced": mode != "off", "passed": passed,
        "reason": None if passed else (f"SA score {sa} is above the '{mode}' limit of {limit}"),
    }


def synth_filter(mol, mode: str = "standard") -> Dict[str, Any]:
    """Gate used by docking (mirrors lipinski_filter): returns the analysis with 'passed'."""
    return analyze(mol, mode)


def synth_rejection_message(v: Dict[str, Any]) -> str:
    return (f"Synthesizability gate: this molecule is rated {v['tier_label'].split(' ', 1)[-1].lower()} to make "
            f"(SA score {v['sa']}, easy = 1, impossible = 10), but '{v['mode']}' mode allows at most {v['limit']}. "
            f"Simplify it (fewer stereocentres, no macrocycles/cages, smaller) or set the Synthesizability "
            f"filter to 'Off' to dock it anyway.")


def synth_summary(smiles_or_mol) -> Optional[Dict[str, Any]]:
    """Light version for lists (similarity cards, guide). None if the molecule cannot be scored."""
    try:
        v = analyze(smiles_or_mol, "off")
    except Exception:
        return None
    if not v.get("available"):
        return {"available": False}
    return {k: v[k] for k in ("available", "sa", "score", "norm", "tier", "tier_label", "colour")} | {"drivers": v["drivers"][:3], "explanation": v["explanation"]}


def composite_score(energy: float, qed: float, synth_norm: Optional[float] = None) -> float:
    """THE game score. Lower energy, higher QED and easier synthesis all score higher."""
    s = -float(energy) * 0.7 + float(qed) * 3.0
    if synth_norm is not None:
        s += SYNTH_WEIGHT * float(synth_norm)
    return round(s, 3)
