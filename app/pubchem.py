# ══ CELL 09 · app/pubchem.py ══
"""
pubchem.py — PubChem lookups with tolerant SMILES-field handling.

PubChem renamed its SMILES properties in 2025
(IsomericSMILES→SMILES, CanonicalSMILES→ConnectivitySMILES). Old code kept
reading the old names, got an empty string, and the game then said
"select a molecule". `pick_smiles` accepts every spelling, and
`fetch_properties` retries with the legacy property names if PubChem answers 400.
"""
import json
import logging
import urllib.error
import urllib.parse
import urllib.request
from typing import Dict, List, Optional

logger = logging.getLogger(__name__)
UA = {"User-Agent": "MedChemGame/3.1"}
BASE = "https://pubchem.ncbi.nlm.nih.gov/rest/pug/compound"

_NEW = "SMILES,ConnectivitySMILES,MolecularFormula,MolecularWeight,IUPACName"
_OLD = "IsomericSMILES,CanonicalSMILES,MolecularFormula,MolecularWeight,IUPACName"


def pick_smiles(p: dict) -> str:
    for key in ("SMILES", "IsomericSMILES", "ConnectivitySMILES", "CanonicalSMILES"):
        v = p.get(key)
        if v:
            return v
    return ""


def _get_json(url: str, timeout: int = 10):
    with urllib.request.urlopen(urllib.request.Request(url, headers=UA), timeout=timeout) as r:
        return json.loads(r.read().decode())


def fetch_properties(namespace: str, identifier: str, timeout: int = 10) -> Optional[Dict]:
    """namespace = 'name' or 'cid'. Returns the first property dict (with 'smiles' added) or None."""
    ident = urllib.parse.quote(str(identifier))
    last_err = None
    for props in (_NEW, _OLD):
        try:
            data = _get_json(f"{BASE}/{namespace}/{ident}/property/{props}/JSON", timeout)
            rows = data.get("PropertyTable", {}).get("Properties", [])
            if rows:
                row = dict(rows[0])
                row["smiles"] = pick_smiles(row)
                return row
            return None
        except urllib.error.HTTPError as e:
            last_err = e
            if e.code == 404:
                return None
            if e.code != 400:
                raise
        except Exception as e:
            last_err = e
            raise
    logger.warning("PubChem property lookup failed for %s=%s: %s", namespace, identifier, last_err)
    return None


def autocomplete(query: str, limit: int = 8) -> List[str]:
    try:
        d = _get_json(f"https://pubchem.ncbi.nlm.nih.gov/rest/autocomplete/compound/{urllib.parse.quote(query)}/json?limit={limit}")
        return d.get("dictionary_terms", {}).get("compound", []) or []
    except Exception as exc:
        logger.warning("PubChem autocomplete failed: %s", exc)
        return []


def result_card(name: str, p: Dict) -> Dict:
    cid = p.get("CID", 0)
    return {
        "name": name,
        "cid": cid,
        "smiles": p.get("smiles", ""),
        "formula": p.get("MolecularFormula", ""),
        "mw": p.get("MolecularWeight", 0),
        "iupac_name": p.get("IUPACName", ""),
        "image_url": f"{BASE}/cid/{cid}/PNG?image_size=200x200" if cid else None,
    }
