"""
oral.py — the drug score and the built-in library of known or marketed drugs.

The drug score (0 to 100) says how drug-like a molecule looks. Tanimoto by itself only measures how ALIKE two
molecules are, so it is used here for one part only: resemblance to known or marketed drugs.

    Lipinski rule of five   molecular weight <= 500, logP <= 5, H-bond donors <= 5, H-bond acceptors <= 10
                            (each rule has its own weight: size and greasiness count most, acceptors least)
    Veber rules             polar surface area <= 140 A^2, rotatable bonds <= 10
    QED                     Bickerton's quantitative estimate of drug-likeness, 0 to 1
    Resemblance             best MACCS Tanimoto to the known or marketed drugs in the library

All of this is a rule-of-thumb screen, not a prediction of how a molecule behaves in a body.
"""
from typing import Any, Dict, List, Optional, Tuple

from rdkit import Chem, DataStructs
from rdkit.Chem import MACCSkeys, QED, rdMolDescriptors

from .lipinski import lipinski_filter

# (name, SMILES, note) — approved drugs taken by mouth
ORAL_DRUGS: List[Tuple[str, str, str]] = [
    ("Aspirin", "CC(=O)Oc1ccccc1C(=O)O", "pain / antiplatelet"),
    ("Paracetamol", "CC(=O)Nc1ccc(O)cc1", "pain / fever"),
    ("Ibuprofen", "CC(C)Cc1ccc(cc1)C(C)C(=O)O", "NSAID"),
    ("Naproxen", "COc1ccc2cc([C@H](C)C(=O)O)ccc2c1", "NSAID"),
    ("Diclofenac", "OC(=O)Cc1ccccc1Nc1c(Cl)cccc1Cl", "NSAID"),
    ("Caffeine", "Cn1cnc2c1c(=O)n(C)c(=O)n2C", "stimulant"),
    ("Metformin", "CN(C)C(=N)NC(N)=N", "type 2 diabetes"),
    ("Propranolol", "CC(C)NCC(O)COc1cccc2ccccc12", "beta blocker"),
    ("Atenolol", "CC(C)NCC(O)COc1ccc(CC(N)=O)cc1", "beta blocker"),
    ("Metoprolol", "COCCc1ccc(OCC(O)CNC(C)C)cc1", "beta blocker"),
    ("Amlodipine", "CCOC(=O)C1=C(COCCN)NC(C)=C(C1c1ccccc1Cl)C(=O)OC", "blood pressure"),
    ("Nifedipine", "COC(=O)C1=C(C)NC(C)=C(C1c1ccccc1[N+](=O)[O-])C(=O)OC", "blood pressure"),
    ("Captopril", "CC(CS)C(=O)N1CCCC1C(O)=O", "ACE inhibitor"),
    ("Lisinopril", "NCCCCC(NC(CCc1ccccc1)C(O)=O)C(=O)N1CCCC1C(O)=O", "ACE inhibitor"),
    ("Losartan", "CCCCc1nc(Cl)c(CO)n1Cc1ccc(cc1)-c1ccccc1-c1nn[nH]n1", "blood pressure"),
    ("Hydrochlorothiazide", "NS(=O)(=O)c1cc2c(NCNS2(=O)=O)cc1Cl", "diuretic"),
    ("Furosemide", "NS(=O)(=O)c1cc(C(O)=O)c(NCc2ccco2)cc1Cl", "diuretic"),
    ("Warfarin", "CC(=O)CC(c1ccccc1)c1c(O)c2ccccc2oc1=O", "anticoagulant"),
    ("Diazepam", "CN1C(=O)CN=C(c2ccccc2)c2cc(Cl)ccc12", "anxiety"),
    ("Fluoxetine", "CNCCC(Oc1ccc(cc1)C(F)(F)F)c1ccccc1", "depression"),
    ("Sertraline", "CNC1CCC(c2ccc(Cl)c(Cl)c2)c2ccccc12", "depression"),
    ("Haloperidol", "OC1(CCN(CCCC(=O)c2ccc(F)cc2)CC1)c1ccc(Cl)cc1", "antipsychotic"),
    ("Risperidone", "Cc1nc2CCCCn2c(=O)c1CCN1CCC(CC1)c1noc2cc(F)ccc12", "antipsychotic"),
    ("Olanzapine", "Cc1cc2C(=Nc3ccccc3Nc2s1)N1CCN(C)CC1", "antipsychotic"),
    ("Carbamazepine", "NC(=O)N1c2ccccc2C=Cc2ccccc12", "epilepsy"),
    ("Phenytoin", "O=C1NC(=O)C(N1)(c1ccccc1)c1ccccc1", "epilepsy"),
    ("Lamotrigine", "Nc1nnc(c(N)n1)-c1cccc(Cl)c1Cl", "epilepsy"),
    ("Gabapentin", "NCC1(CC(O)=O)CCCCC1", "nerve pain"),
    ("Tramadol", "COc1cccc(c1)C1(O)CCCCC1CN(C)C", "pain"),
    ("Codeine", "COc1ccc2CC3N(C)CCC45C(Oc1c24)C(O)C=CC35", "pain"),
    ("Morphine", "CN1CCC23C4C1CC5=C2C(=C(C=C5)O)OC3C(C=C4)O", "pain"),
    ("Omeprazole", "COc1ccc2[nH]c(nc2c1)S(=O)Cc1ncc(C)c(OC)c1C", "acid reflux"),
    ("Ranitidine", "CNC(NCCSCc1ccc(CN(C)C)o1)=C[N+]([O-])=O", "acid reflux"),
    ("Ciprofloxacin", "OC(=O)C1=CN(C2CC2)c2cc(N3CCNCC3)c(F)cc2C1=O", "antibiotic"),
    ("Amoxicillin", "CC1(C)SC2C(NC(=O)C(N)c3ccc(O)cc3)C(=O)N2C1C(O)=O", "antibiotic"),
    ("Sildenafil", "CCCc1nn(C)c2c1nc([nH]c2=O)-c1cc(ccc1OCC)S(=O)(=O)N1CCN(C)CC1", "erectile dysfunction"),
    ("Loratadine", "CCOC(=O)N1CCC(=C2c3ccc(Cl)cc3CCc3cccnc23)CC1", "allergy"),
    ("Cetirizine", "OC(=O)COCCN1CCN(CC1)C(c1ccccc1)c1ccc(Cl)cc1", "allergy"),
    ("Ondansetron", "Cc1nccn1CC1CCc2c(C1=O)c1ccccc1n2C", "nausea"),
    ("Prednisone", "CC12CC(=O)C3C(CCC4=CC(=O)C=CC34C)C1CCC2(O)C(=O)CO", "steroid"),
    ("Sumatriptan", "CNS(=O)(=O)Cc1ccc2[nH]cc(CCN(C)C)c2c1", "migraine"),
    ("Zolpidem", "Cc1ccc(cc1)-c1nc2ccc(C)cn2c1CC(=O)N(C)C", "insomnia"),
    ("Selegiline", "CC(Cc1ccccc1)N(C)CC#C", "Parkinson's (MAO-B)"),
    ("Rasagiline", "C#CCNC1CCc2ccccc21", "Parkinson's (MAO-B)"),
    ("Safinamide", "CC(NCc1ccc(OCc2cccc(F)c2)cc1)C(N)=O", "Parkinson's (MAO-B)"),
    ("Moclobemide", "O=C(NCCN1CCOCC1)c1ccc(Cl)cc1", "depression (MAO-A)"),
    ("Ketoconazole", "CC(=O)N1CCN(CC1)c1ccc(OCC2COC(Cn3ccnc3)(O2)c2ccc(Cl)cc2Cl)cc1", "antifungal (CYP3A4 inhibitor)"),
    ("Midazolam", "Cc1ncc2CN=C(c3ccccc3F)c3cc(Cl)ccc3-n12", "sedative (CYP3A4 substrate)"),
    ("Donepezil", "COc1cc2CC(CC3CCN(Cc4ccccc4)CC3)C(=O)c2cc1OC", "Alzheimer's"),
    ("Simvastatin", "CCC(C)(C)C(=O)OC1CC(C)C=C2C=CC(C)C(CCC3CC(O)CC(=O)O3)C12", "cholesterol"),
    ("Glibenclamide", "COc1ccc(Cl)cc1C(=O)NCCc1ccc(cc1)S(=O)(=O)NC(=O)NC1CCCCC1", "type 2 diabetes"),
    ("Pioglitazone", "CCc1ccc(CCOc2ccc(CC3SC(=O)NC3=O)cc2)nc1", "type 2 diabetes"),
]

# Known ligands per built-in target: used as the yardstick ("reference docking") the player has to beat.
TARGET_LIGANDS: Dict[str, List[Tuple[str, str, str]]] = {
    "MAOB": [("Selegiline", "CC(Cc1ccccc1)N(C)CC#C", "approved MAO-B inhibitor"),
             ("Rasagiline", "C#CCNC1CCc2ccccc21", "approved MAO-B inhibitor"),
             ("Safinamide", "CC(NCc1ccc(OCc2cccc(F)c2)cc1)C(N)=O", "approved MAO-B inhibitor")],
    "MAOA": [("Clorgyline", "C#CCN(C)CCCOc1ccc(Cl)cc1Cl", "classic MAO-A inhibitor (research)"),
             ("Moclobemide", "O=C(NCCN1CCOCC1)c1ccc(Cl)cc1", "approved MAO-A inhibitor"),
             ("Harmine", "COc1ccc2c(c1)[nH]c1c(C)nccc12", "natural MAO-A inhibitor")],
    "CYP3A4": [("Ketoconazole", "CC(=O)N1CCN(CC1)c1ccc(OCC2COC(Cn3ccnc3)(O2)c2ccc(Cl)cc2Cl)cc1", "azole CYP3A4 inhibitor"),
               ("Midazolam", "Cc1ncc2CN=C(c3ccccc3F)c3cc(Cl)ccc3-n12", "classic CYP3A4 substrate"),
               ("Nifedipine", "COC(=O)C1=C(C)NC(C)=C(C1c1ccccc1[N+](=O)[O-])C(=O)OC", "CYP3A4 substrate"),
               ("Testosterone", "CC12CCC3C(CCC4=CC(=O)CCC34C)C1CCC2O", "CYP3A4 substrate")],
    "BACE1": [("Verubecestat", "CC1(CS(=O)(=O)N(C)C(N)=N1)c1cc(NC(=O)c2ccc(F)cn2)ccc1F", "clinical BACE1 inhibitor (amidine head group)")],
}

_cache: Dict[str, Any] = {}


def _fp(m):
    return MACCSkeys.GenMACCSKeys(m)


def _library() -> List[Dict[str, Any]]:
    """Fingerprinted reference library (oral drugs + target ligands), built once."""
    if "lib" not in _cache:
        rows, seen = [], set()
        def add(name, smi, note, kind, target=None):
            m = Chem.MolFromSmiles(smi)
            if m is None:
                return
            canon = Chem.MolToSmiles(m)
            for r in rows:
                if r["canon"] == canon:
                    if target and target not in r["targets"]:
                        r["targets"].append(target)
                    if kind == "oral":
                        r["oral"] = True
                    return
            rows.append({"name": name, "smiles": smi, "note": note, "canon": canon, "fp": _fp(m),
                         "oral": kind == "oral", "targets": [target] if target else []})
        for n, s, note in ORAL_DRUGS:
            add(n, s, note, "oral")
        for tid, lst in TARGET_LIGANDS.items():
            for n, s, note in lst:
                add(n, s, note, "target", tid)
        _cache["lib"] = rows
    return _cache["lib"]


def reference_library() -> List[Tuple[str, str]]:
    """(name, smiles) pairs for the Drug Twin Finder's always-available 'known drugs' source."""
    return [(r["name"], r["smiles"]) for r in _library()]


def nearest_known(mol, top: int = 5, oral_only: bool = False) -> List[Dict[str, Any]]:
    q = _fp(mol)
    out = []
    for r in _library():
        if oral_only and not r["oral"]:
            continue
        t = DataStructs.TanimotoSimilarity(q, r["fp"])
        out.append({"name": r["name"], "smiles": r["smiles"], "note": r["note"], "maccs": round(t, 4),
                    "percent": round(t * 100), "targets": r["targets"], "oral": r["oral"]})
    out.sort(key=lambda d: d["maccs"], reverse=True)
    return out[:top]


# ─────────────────────────── the drug score ───────────────────────────
# Weights of the four parts (must add up to 100). Resemblance gets the least because the library has only 56 molecules.
DRUG_WEIGHTS = {"ro5": 35, "veber": 20, "qed": 30, "sim": 15}
# Lipinski rules: (key, label, limit, weight inside the Lipinski part, unit). The weights add up to 100.
#   molecular weight and logP decide both dissolving and crossing membranes -> highest weight
#   H-bond donors cost the most energy to strip of water when entering a membrane -> next
#   H-bond acceptors cost the least -> lowest weight
LIPINSKI_RULES = [("mw", "Molecular weight", 500.0, 30, "Da"), ("logp", "logP (greasiness)", 5.0, 30, ""),
                  ("hbd", "H-bond donors", 5, 25, ""), ("hba", "H-bond acceptors", 10, 15, "")]
VEBER_RULES = [("tpsa", "Polar surface area", 140.0, 55, "Å²"), ("rot", "Rotatable bonds", 10, 45, "")]
OVERSHOOT = 0.5      # a rule loses its points gradually past its limit and reaches zero credit 50% beyond it


def _credit(value: float, limit: float) -> float:
    if value <= limit:
        return 1.0
    return max(0.0, 1.0 - (value - limit) / (limit * OVERSHOOT))


def _weighted(values: Dict[str, float], rules) -> Tuple[int, List[Dict[str, Any]]]:
    pts, sub = 0.0, []
    for key, label, limit, weight, unit in rules:
        v = values[key]
        c = _credit(v, limit)
        pts += weight * c
        sub.append({"key": key, "label": label, "value": round(float(v), 2), "limit": limit, "unit": unit,
                    "weight": weight, "credit": round(c, 2), "ok": v <= limit})
    return int(round(pts)), sub


def _grade(score: int) -> Tuple[str, str]:
    if score >= 80: return "excellent", "🟢 Excellent drug-like profile"
    if score >= 60: return "good", "🟢 Good drug-like profile"
    if score >= 40: return "borderline", "🟡 Borderline"
    return "poor", "🔴 Poor drug-likeness"


def _over(sub) -> str:
    bad = [s for s in sub if not s["ok"]]
    if not bad:
        return ""
    return "; ".join(f"{s['label']} {s['value']:g} is over {s['limit']:g}{(' ' + s['unit']) if s['unit'] else ''} "
                     f"(keeps {round(s['credit'] * 100)}% of its {s['weight']} points)" for s in bad)


def drug_profile(mol: Chem.Mol) -> Dict[str, Any]:
    """The drug score: how drug-like the molecule looks, from 0 to 100, with a transparent breakdown."""
    m = Chem.RemoveHs(mol)
    ro5 = lipinski_filter(m, mode="standard")
    pr = ro5["properties"]
    ro5_pts, ro5_sub = _weighted({"mw": pr["mw"], "logp": pr["logp"], "hbd": pr["hbd"], "hba": pr["hba"]}, LIPINSKI_RULES)
    rot = int(rdMolDescriptors.CalcNumRotatableBonds(m))
    tpsa = float(rdMolDescriptors.CalcTPSA(m))
    veb_pts, veb_sub = _weighted({"tpsa": tpsa, "rot": rot}, VEBER_RULES)
    try:
        qed = float(QED.qed(m))
    except Exception:
        qed = 0.0
    near = nearest_known(m, top=3)                       # all known or marketed drugs in the library
    best = near[0] if near else {"maccs": 0.0, "name": "—", "percent": 0}
    sim_pts = round(min(1.0, best["maccs"] / 0.8) * 100)
    parts = [
        {"key": "ro5", "label": "Lipinski rule of five", "weight": DRUG_WEIGHTS["ro5"], "points": ro5_pts, "sub": ro5_sub,
         "detail": _over(ro5_sub) or "all four limits met", "ok": ro5["violations"] <= 1},
        {"key": "veber", "label": "Veber rules (flexibility, polarity)", "weight": DRUG_WEIGHTS["veber"], "points": veb_pts, "sub": veb_sub,
         "detail": _over(veb_sub) or f"{rot} rotatable bonds (limit 10), polar surface area {tpsa:.0f} Å² (limit 140)",
         "ok": all(s["ok"] for s in veb_sub)},
        {"key": "qed", "label": "QED drug-likeness", "weight": DRUG_WEIGHTS["qed"], "points": round(qed * 100), "sub": [],
         "detail": f"QED {qed:.2f} (0 = unlike drugs, 1 = very drug-like)", "ok": qed >= 0.5},
        {"key": "sim", "label": "Resemblance to known or marketed drugs", "weight": DRUG_WEIGHTS["sim"], "points": sim_pts, "sub": [],
         "detail": f"closest: {best['name']} ({best['percent']}% alike by MACCS Tanimoto; 80% or more earns full marks)",
         "ok": best["maccs"] >= 0.5},
    ]
    score = int(round(sum(p["weight"] * p["points"] for p in parts) / 100))
    grade, grade_label = _grade(score)
    return {"score": score, "grade": grade, "grade_label": grade_label, "parts": parts, "weights": dict(DRUG_WEIGHTS),
            "lipinski": ro5, "veber": {"passed": all(s["ok"] for s in veb_sub), "rot_bonds": rot, "tpsa": round(tpsa, 1)},
            "qed": round(qed, 3), "nearest": near, "nearest_oral": near,
            "note": ("A rule-of-thumb screen, not a prediction of how a molecule behaves in a body. Each part is scored 0 to 100 "
                     "and mixed with the weights shown. Resemblance compares with the 56 known or marketed drugs in the library; "
                     "because that list is short, it carries the smallest weight.")}


oral_profile = drug_profile            # old name kept so nothing that imports it breaks


# ─────────────────────────── comparison molecules for a (protein, ligand) pair ───────────────────────────
MAX_REFERENCES = 3


def plan_references(mol, target_id: Optional[str]) -> Dict[str, Any]:
    """Which molecules to dock against this protein to build the score to beat for THIS ligand:
         * the ligand itself, if it is a known molecule (its own best known score for this protein),
         * the molecules in the library most similar to the ligand,
         * the known binder of this protein most similar to the ligand."""
    q = _fp(mol)
    mycanon = Chem.MolToSmiles(mol)
    scored = []
    for r in _library():
        t = DataStructs.TanimotoSimilarity(q, r["fp"])
        scored.append({"name": r["name"], "smiles": r["smiles"], "canon": r["canon"], "note": r["note"], "maccs": round(t, 4),
                       "percent": round(t * 100), "targets": list(r["targets"]),
                       "binder": bool(target_id) and target_id in r["targets"]})
    scored.sort(key=lambda d: d["maccs"], reverse=True)
    same = next((s for s in scored if s["canon"] == mycanon), None)
    binder = next((s for s in scored if s["binder"] and s["canon"] != mycanon), None)
    chosen = [dict(same, kind="same-molecule")] if same else []
    room = MAX_REFERENCES - len(chosen) - (1 if binder else 0)          # keep one place for the known binder
    for s in scored:
        if room <= 0:
            break
        if s["canon"] == mycanon or (binder and s["canon"] == binder["canon"]):
            continue
        chosen.append(dict(s, kind="similar-drug"))
        room -= 1
    if binder:
        chosen.append(dict(binder, kind="known-binder"))
    chosen.sort(key=lambda c: (c["kind"] != "same-molecule", -c["maccs"]))
    others = [c for c in chosen if c["kind"] != "same-molecule"]
    top = max([c["maccs"] for c in others], default=0.0)
    loose = (not same) and top < 0.55
    if same:
        why = (f"{same['name']} is a known molecule, so its own docking score with this protein is the record. "
               "Modify it and try to beat that score.")
    elif loose:
        why = ("Nothing in the library is very close to your molecule (best match "
               f"{round(top * 100)}%), so the score to beat is only a loose yardstick.")
    else:
        why = ("The score to beat is the best docking score, for this protein, among molecules similar to yours "
               f"(closest match {round(top * 100)}%)" + (" and the closest known binder of this protein." if binder else "."))
    return {"target_id": target_id, "candidates": chosen, "loose": loose, "same": bool(same), "rationale": why}


def pick_reference(mol, target_id: Optional[str]) -> Dict[str, Any]:
    """Older single-reference view of the plan (kept for compatibility)."""
    plan = plan_references(mol, target_id)
    first = plan["candidates"][0]
    return {"source": "target-known" if first["kind"] == "known-binder" else "nearest-drug", "reference": first,
            "alternatives": plan["candidates"], "rationale": plan["rationale"]}



# BEGIN LIVE_EXTERNAL_REFERENCE_PATCH
import json as _live_json
import sys as _live_sys
import time as _live_time
import urllib.parse as _live_parse
import urllib.request as _live_request

_LIVE_REFERENCE_CACHE = {}
_LIVE_HEADERS = {"User-Agent": "MedChemGame/3.1", "Accept": "application/json"}


def _live_sa_score(mol):
    try:
        from .synth import analyze as _synth_analyze
        report = _synth_analyze(mol, "off")
        return report.get("sa") if isinstance(report, dict) else None
    except Exception:
        return None


def _live_similar(smiles, min_similarity=0.35, max_records=30):
    """Live ChEMBL + PubChem similarity lookup. Never uses ORAL_DRUGS."""
    mol = Chem.MolFromSmiles((smiles or "").strip())
    if mol is None:
        return []

    canonical = Chem.MolToSmiles(mol)
    key = (canonical, round(float(min_similarity), 2), int(max_records))

    cached = _LIVE_REFERENCE_CACHE.get(key)
    if cached and _live_time.time() - cached["time"] < 1800:
        return [dict(row) for row in cached["rows"]]

    query_fp = _fp(mol)
    encoded = _live_parse.quote(canonical, safe="")
    limit = max(1, min(int(max_records), 50))
    found = {}

    def add(name, smi, source, note="", chembl="", cid=""):
        hit = Chem.MolFromSmiles(smi or "")
        if hit is None:
            return
        hit_canonical = Chem.MolToSmiles(hit)
        if hit_canonical == canonical:
            return
        found.setdefault(hit_canonical, {
            "name": name or (f"PubChem CID {cid}" if cid else "Database compound"),
            "smiles": smi,
            "canon": hit_canonical,
            "note": note or source,
            "source": source,
            "chembl": chembl,
            "cid": cid,
        })

    # ChEMBL live similarity API
    try:
        threshold = max(40, min(95, int(float(min_similarity) * 100)))
        url = (
            "https://www.ebi.ac.uk/chembl/api/data/similarity/"
            f"{encoded}/{threshold}.json?limit={limit}"
        )
        request = _live_request.Request(url, headers=_LIVE_HEADERS)
        with _live_request.urlopen(request, timeout=15) as response:
            data = _live_json.loads(response.read().decode("utf-8"))

        for item in data.get("molecules", []) or []:
            smi = (item.get("molecule_structures") or {}).get("canonical_smiles")
            if smi:
                phase = item.get("max_phase")
                add(
                    item.get("pref_name") or item.get("molecule_chembl_id"),
                    smi, "chembl",
                    f"ChEMBL{f' · Phase {phase}' if phase not in (None, '') else ''}",
                    item.get("molecule_chembl_id", "")
                )
    except Exception:
        pass

    # PubChem live similarity API — valid fastsimilarity_2d SMILES endpoint
    try:
        threshold = max(50, min(95, int(float(min_similarity) * 100)))
        url = (
            "https://pubchem.ncbi.nlm.nih.gov/rest/pug/compound/"
            f"fastsimilarity_2d/smiles/{encoded}/cids/JSON"
            f"?Threshold={threshold}&MaxRecords={limit}"
        )
        request = _live_request.Request(url, headers=_LIVE_HEADERS)
        with _live_request.urlopen(request, timeout=20) as response:
            data = _live_json.loads(response.read().decode("utf-8"))

        cids = (data.get("IdentifierList") or {}).get("CID", [])
        if cids:
            cid_text = ",".join(str(cid) for cid in cids[:limit])
            prop_url = (
                "https://pubchem.ncbi.nlm.nih.gov/rest/pug/compound/cid/"
                f"{cid_text}/property/Title,SMILES,ConnectivitySMILES/JSON"
            )
            request = _live_request.Request(prop_url, headers=_LIVE_HEADERS)
            with _live_request.urlopen(request, timeout=20) as response:
                properties = _live_json.loads(response.read().decode("utf-8"))

            for item in properties.get("PropertyTable", {}).get("Properties", []):
                smi = item.get("SMILES") or item.get("ConnectivitySMILES")
                if smi:
                    cid = str(item.get("CID") or "")
                    add(item.get("Title") or f"PubChem CID {cid}", smi, "pubchem", f"PubChem CID {cid}", "", cid)
    except Exception:
        pass

    rows = []
    for row in found.values():
        hit = Chem.MolFromSmiles(row["smiles"])
        if hit is None:
            continue
        similarity = float(DataStructs.TanimotoSimilarity(query_fp, _fp(hit)))
        rows.append({
            **row,
            "maccs": round(similarity, 4),
            "percent": round(similarity * 100),
            "targets": [],
            "oral": False,
        })

    rows.sort(key=lambda row: row["maccs"], reverse=True)
    rows = rows[:limit]

    for row in rows:
        hit = Chem.MolFromSmiles(row["smiles"])
        row["sa_score"] = _live_sa_score(hit) if hit is not None else None

    _LIVE_REFERENCE_CACHE[key] = {"time": _live_time.time(), "rows": rows}
    return [dict(row) for row in rows]


_old_nearest_known = nearest_known
_old_drug_profile = drug_profile
_old_plan_references = plan_references
_old_reference_library = reference_library


def nearest_known(mol, top=5, oral_only=False):
    try:
        return _live_similar(Chem.MolToSmiles(mol), min_similarity=0.35, max_records=max(20, int(top) * 5))[:top]
    except Exception:
        return []


def drug_profile(mol):
    # The original Lipinski, Veber, QED, and score math stays unchanged.
    # Its resemblance component now calls the live nearest_known above.
    result = _old_drug_profile(mol)
    nearest = nearest_known(mol, top=3)
    result["nearest"] = nearest
    result["nearest_oral"] = nearest
    result["note"] = (
        "Drug score combines Lipinski, Veber and QED with resemblance to live "
        "ChEMBL and PubChem molecules. The 52 built-in Cell 49o molecules are not used."
    )

    for part in result.get("parts", []):
        if part.get("key") == "sim":
            best = nearest[0] if nearest else None
            part["label"] = "Resemblance to ChEMBL / PubChem molecules"
            part["detail"] = (
                f"Closest live result: {best['name']} ({best['percent']}% MACCS similarity)"
                if best else "No live ChEMBL or PubChem reference was available."
            )
            part["ok"] = bool(best and best["maccs"] >= 0.5)
    return result


def plan_references(mol, target_id=None):
    """Reference docking candidates from ChEMBL/PubChem only."""
    rows = _live_similar(Chem.MolToSmiles(mol), min_similarity=0.35, max_records=30)
    maximum = max(1, int(globals().get("MAX_REFERENCES", 3)))

    candidates = []
    for row in rows[:maximum]:
        candidates.append({
            **row,
            "binder": False,
            "kind": "similar-drug",
            "canon": row["canon"],
        })

    best_percent = candidates[0]["percent"] if candidates else 0
    return {
        "target_id": target_id,
        "candidates": candidates,
        "same": False,
        "loose": (not candidates) or best_percent < 55,
        "rationale": (
            f"Reference docking uses the closest live ChEMBL/PubChem molecules "
            f"(best similarity: {best_percent}%)."
            if candidates else
            "No live ChEMBL or PubChem reference is currently available."
        ),
    }


def reference_library():
    # Prevent Drug Twin Finder from adding the 52 Cell 49o molecules.
    return []


def _rebind_imports(old, new):
    this_module = _live_sys.modules.get(__name__)
    for module in list(_live_sys.modules.values()):
        if module is None or module is this_module or not getattr(module, "__name__", "").startswith("app"):
            continue
        for name, value in list(getattr(module, "__dict__", {}).items()):
            if value is old:
                setattr(module, name, new)


_rebind_imports(_old_nearest_known, nearest_known)
_rebind_imports(_old_drug_profile, drug_profile)
_rebind_imports(_old_plan_references, plan_references)
_rebind_imports(_old_reference_library, reference_library)
oral_profile = drug_profile
# END LIVE_EXTERNAL_REFERENCE_PATCH



# BEGIN TWO_REFERENCE_PLAN_PATCH
MAX_REFERENCES = 2
_two_ref_old_plan = plan_references

def _closest_target_binder(mol, target_id):
    choices = TARGET_LIGANDS.get((target_id or "").upper(), [])
    query_fp = _fp(mol)
    best = None

    for name, smiles, note in choices:
        binder_mol = Chem.MolFromSmiles(smiles)
        if binder_mol is None:
            continue
        similarity = float(DataStructs.TanimotoSimilarity(query_fp, _fp(binder_mol)))
        candidate = {
            "name": name,
            "smiles": smiles,
            "canon": Chem.MolToSmiles(binder_mol),
            "note": note,
            "source": "known-target-binder",
            "maccs": round(similarity, 4),
            "percent": round(similarity * 100),
            "targets": [(target_id or "").upper()],
            "binder": True,
            "kind": "known-binder",
            "sa_score": _live_sa_score(binder_mol) if "_live_sa_score" in globals() else None,
        }
        if best is None or candidate["maccs"] > best["maccs"]:
            best = candidate
    return best


def plan_references(mol, target_id=None):
    # One target-specific closest binder + one live external database analogue.
    binder = _closest_target_binder(mol, target_id)
    live = _live_similar(Chem.MolToSmiles(mol), min_similarity=0.35, max_records=20)

    candidates = []
    if binder:
        candidates.append(binder)

    used = {binder["canon"]} if binder else set()
    for row in live:
        if row["canon"] in used:
            continue
        candidates.append({
            **row,
            "binder": False,
            "kind": "similar-drug",
        })
        break

    candidates = candidates[:2]

    return {
        "target_id": target_id,
        "candidates": candidates,
        "same": False,
        "loose": not bool(candidates),
        "rationale": (
            "Two docking references only: the closest known binder for this protein "
            "and the closest live ChEMBL/PubChem analogue. Their scores use the same engine and grid."
            if candidates else
            "No live reference compound is currently available."
        ),
    }


# Update modules that imported the older plan function.
import sys as _two_ref_sys
for _module in list(_two_ref_sys.modules.values()):
    if _module is None or getattr(_module, "__name__", None) == __name__:
        continue
    if getattr(_module, "__name__", "").startswith("app"):
        for _name, _value in list(getattr(_module, "__dict__", {}).items()):
            if _value is _two_ref_old_plan:
                setattr(_module, _name, plan_references)
# END TWO_REFERENCE_PLAN_PATCH

# BEGIN UNIQUE_REFERENCE_PLAN_V1
MAX_REFERENCES = 2
_previous_reference_plan_v1 = plan_references

def _distinct_reference_plan_v1(mol, target_id=None):
    query_smiles = Chem.MolToSmiles(mol, isomericSmiles=True)
    query_fp = _fp(mol)
    candidates = []
    selected_fps = []
    selected_canons = set()

    def add(candidate, kind):
        if not candidate:
            return False

        smiles = candidate.get("smiles", "")
        hit = Chem.MolFromSmiles(smiles)
        if hit is None:
            return False

        canonical = Chem.MolToSmiles(hit, isomericSmiles=True)
        fp = _fp(hit)
        similarity_to_query = float(DataStructs.TanimotoSimilarity(query_fp, fp))

        # Do not dock the loaded molecule, a renamed salt, or a 100%-identical twin.
        if canonical == query_smiles or similarity_to_query >= 0.999999:
            return False

        if canonical in selected_canons:
            return False

        # Do not choose two references which are 100% alike to each other.
        if any(
            DataStructs.TanimotoSimilarity(fp, old_fp) >= 0.999999
            for old_fp in selected_fps
        ):
            return False

        candidates.append({
            **candidate,
            "canon": canonical,
            "kind": kind,
            "binder": kind == "known-binder",
            "maccs": round(similarity_to_query, 4),
            "percent": round(similarity_to_query * 100),
        })
        selected_canons.add(canonical)
        selected_fps.append(fp)
        return True

    # Preserve the existing target-specific binder when it is not identical.
    try:
        add(_closest_target_binder(mol, target_id), "known-binder")
    except Exception:
        pass

    # Then choose genuinely distinct database analogues until there are two.
    try:
        live_rows = _live_similar(query_smiles, min_similarity=0.35, max_records=50)
    except Exception:
        live_rows = []

    for row in live_rows:
        add(row, "similar-drug")
        if len(candidates) >= MAX_REFERENCES:
            break

    return {
        "target_id": target_id,
        "candidates": candidates,
        "same": False,
        "loose": not bool(candidates),
        "rationale": (
            "Up to two distinct comparison molecules are used. "
            "Exact and 100%-similar duplicates are skipped automatically."
            if candidates else
            "No sufficiently distinct comparison molecule is currently available."
        ),
    }


def plan_references(mol, target_id=None):
    return _distinct_reference_plan_v1(mol, target_id)


# Update modules which imported plan_references before this patch.
import sys as _unique_ref_sys
for _module in list(_unique_ref_sys.modules.values()):
    if _module is None or not getattr(_module, "__name__", "").startswith("app"):
        continue
    for _name, _value in list(getattr(_module, "__dict__", {}).items()):
        if _value is _previous_reference_plan_v1:
            setattr(_module, _name, plan_references)
# END UNIQUE_REFERENCE_PLAN_V1


# BEGIN DISTINCT_REFERENCES_V2
MAX_REFERENCES = 2
_previous_reference_plan_v2 = plan_references

def plan_references(mol, target_id=None):
    query_smiles = Chem.MolToSmiles(mol, isomericSmiles=True)
    query_fp = _fp(mol)
    chosen = []
    chosen_fps = []
    chosen_canons = set()

    def add(candidate, kind):
        if not candidate:
            return False

        candidate_mol = Chem.MolFromSmiles(candidate.get("smiles", ""))
        if candidate_mol is None:
            return False

        canonical = Chem.MolToSmiles(candidate_mol, isomericSmiles=True)
        fp = _fp(candidate_mol)
        sim_to_query = float(DataStructs.TanimotoSimilarity(query_fp, fp))

        # Do not reference-dock the loaded molecule or a 100%-similar duplicate.
        if canonical == query_smiles or sim_to_query >= 0.999999:
            return False

        # Do not use a second reference that is 100% alike to one already chosen.
        if canonical in chosen_canons or any(
            DataStructs.TanimotoSimilarity(fp, previous_fp) >= 0.999999
            for previous_fp in chosen_fps
        ):
            return False

        chosen.append({
            **candidate,
            "canon": canonical,
            "kind": kind,
            "binder": kind == "known-binder",
            "maccs": round(sim_to_query, 4),
            "percent": round(sim_to_query * 100),
        })
        chosen_fps.append(fp)
        chosen_canons.add(canonical)
        return True

    # Existing target-specific binder, only if it is genuinely different.
    try:
        add(_closest_target_binder(mol, target_id), "known-binder")
    except Exception:
        pass

    # Choose the next truly distinct live database analogue.
    try:
        live = _live_similar(query_smiles, min_similarity=0.35, max_records=50)
    except Exception:
        live = []

    for candidate in live:
        add(candidate, "similar-drug")
        if len(chosen) >= MAX_REFERENCES:
            break

    return {
        "target_id": target_id,
        "candidates": chosen,
        "same": False,
        "loose": not bool(chosen),
        "rationale": (
            "Only distinct comparison molecules are reference-docked; "
            "100%-similar duplicates are skipped."
            if chosen else
            "No sufficiently distinct comparison molecule is currently available."
        ),
    }


import sys as _distinct_ref_sys
for _module in list(_distinct_ref_sys.modules.values()):
    if _module is None or not getattr(_module, "__name__", "").startswith("app"):
        continue
    for _name, _value in list(getattr(_module, "__dict__", {}).items()):
        if _value is _previous_reference_plan_v2:
            setattr(_module, _name, plan_references)
# END DISTINCT_REFERENCES_V2
