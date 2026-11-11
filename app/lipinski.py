"""
lipinski.py — Lipinski's Rule of Five (oral drug-likeness rule of thumb).

    molecular weight  <= 500 Da      logP <= 5      H-bond donors <= 5      H-bond acceptors <= 10

Modes:  standard = at most 1 miss allowed   |   strict = no misses allowed   |   off = report only, never reject
This is a screening heuristic: many successful drugs break it, and passing it does not make a molecule a drug.
"""
from typing import Any, Dict, Optional

from rdkit import Chem
from rdkit.Chem import Descriptors

LIPINSKI_LIMITS = {"mw": 500.0, "logp": 5.0, "hbd": 5, "hba": 10}
LIPINSKI_MODES = ("standard", "strict", "off")


def lipinski_filter(mol: Chem.Mol, mode: str = "standard", max_violations: Optional[int] = None) -> Dict[str, Any]:
    if mol is None:
        raise ValueError("Cannot apply the Lipinski filter to an empty molecule")
    mode = (mode or "standard").strip().lower()
    if mode not in LIPINSKI_MODES:
        raise ValueError("Lipinski mode must be 'standard', 'strict' or 'off'")
    if max_violations is None:
        max_violations = {"strict": 0, "standard": 1, "off": 4}[mode]
    max_violations = max(0, int(max_violations))

    m = Chem.RemoveHs(mol)
    values = {"mw": float(Descriptors.MolWt(m)), "logp": float(Descriptors.MolLogP(m)),
              "hbd": int(Descriptors.NumHDonors(m)), "hba": int(Descriptors.NumHAcceptors(m))}
    checks = {k: values[k] <= LIPINSKI_LIMITS[k] for k in ("mw", "logp", "hbd", "hba")}
    messages = {
        "mw": f"molecular weight {values['mw']:.0f} Da > 500 Da",
        "logp": f"logP {values['logp']:.2f} > 5",
        "hbd": f"{values['hbd']} H-bond donors > 5",
        "hba": f"{values['hba']} H-bond acceptors > 10",
    }
    detail = [messages[k] for k in ("mw", "logp", "hbd", "hba") if not checks[k]]
    return {
        "passed": len(detail) <= max_violations,
        "enforced": mode != "off",
        "mode": mode,
        "max_violations": max_violations,
        "violations": len(detail),
        "violations_detail": detail,
        "checks": checks,
        "properties": {"mw": round(values["mw"], 2) + 0.0, "logp": round(values["logp"], 2) + 0.0,
                       "hbd": values["hbd"], "hba": values["hba"]},
    }


def rejection_message(verdict: Dict[str, Any]) -> str:
    n = verdict["violations"]
    return (f"Lipinski gate: this molecule misses {n} Rule-of-Five limit{'s' if n != 1 else ''} "
            f"({'; '.join(verdict['violations_detail'])}), but '{verdict['mode']}' mode allows at most "
            f"{verdict['max_violations']}. Set the Lipinski filter to 'Off' to dock it anyway.")
