
"""
drug_warnings.py — the drug warning system, for any kind of drug and any place it has to act.

No score. RDKit descriptors are compared with thresholds, and the thresholds depend on WHERE the drug must act:
a TPSA of 120 is acceptable for a gut-absorbed drug but is a warning for a brain drug. Every warning is one sentence
(property, value, limit, consequence). Limits are medicinal-chemistry rules of thumb, not hard laws.
"""
from typing import Any, Dict, List, Optional

from rdkit import Chem
from rdkit.Chem import Descriptors

SITES: Dict[str, str] = {
    "oral": "Oral systemic (swallowed; must reach blood)",
    "cns": "Brain / CNS (must cross the blood-brain barrier)",
    "iv": "Intravenous systemic (injected into blood)",
    "topical": "Topical skin (local action)",
    "transdermal": "Transdermal (cross skin into blood)",
    "inhaled": "Inhaled lung (local action)",
    "inhaled_systemic": "Inhaled systemic (lung-to-blood)",
    "ocular": "Ocular (eye drops / local eye action)",
    "intrathecal": "Intrathecal (delivered into cerebrospinal fluid)",
    "gut_local": "Gut-restricted (act in intestine; minimal blood exposure)",
    "buccal": "Buccal / sublingual (mouth lining into blood)",
    "rectal_vaginal": "Rectal / vaginal (local or systemic)",
}
SITE_ORDER = list(SITES)

# Each rule: rdkit name, friendly name, direction ("above" warns when value > limit, "below" when value < limit),
#   limits per site ("*" = every site), unit, group, severe sites (where the warning blocks use), consequence text.
R = {}
def rule(rd, name, direction, limits, unit, group, consequence, severe=()):
    R.setdefault("list", []).append((rd, name, direction, limits, unit, group, consequence, set(severe)))

GROUP = {"size": "Size", "shape": "Shape", "flexibility": "Flexibility", "lipophilicity": "Lipophilicity and solubility",
         "polarity": "Polarity and permeability", "charge": "Charge", "ionisation": "Ionisation",
         "alert": "Structural alerts (toxicity, reactivity)", "metabolism": "Metabolism and clearance",
         "complexity": "Synthesis complexity", "ring": "Ring composition", "surface": "Surface area",
         "electronic": "Electronic", "drug-likeness": "Drug-likeness", "fragment": "Repeated structural features"}

# ── size ──
rule("MolWt", "Molecular weight", "above", {"oral": 500, "cns": 450, "iv": 800, "topical": 500, "transdermal": 500, "inhaled": 400, "inhaled_systemic": 500, "ocular": 500, "intrathecal": 800, "gut_local": 800, "buccal": 500, "rectal_vaginal": 500}, "Da", "size",
     {"default": "Heavy molecules are harder to absorb and to move through tissue.", "oral": "Too heavy to be absorbed well from the gut.",
      "cns": "Too heavy to diffuse into the brain efficiently.", "iv": "Very large: may need a special delivery or formulation.",
      "inhaled": "Too heavy to reach deep into the lungs."}, severe={"oral", "cns"})
rule("HeavyAtomCount", "Heavy atom count", "above", {"oral": 38, "cns": 32, "iv": 60, "topical": 38, "transdermal": 38, "inhaled": 30, "inhaled_systemic": 38, "ocular": 38, "intrathecal": 60, "gut_local": 50, "buccal": 38, "rectal_vaginal": 45}, "atoms", "size",
     "A very large molecule: harder to make, absorb and distribute.")
rule("HeavyAtomCount", "Heavy atom count (too small)", "below", {"*": 12}, "atoms", "size",
     "Very small molecules usually bind too weakly to act as drugs.")
rule("RingCount", "Ring count", "above", {"*": 6}, "", "size", "Many rings often means a hard-to-make, poorly soluble molecule.")
rule("Kappa1", "Shape index Kappa1", "above", {"*": 9}, "", "shape", "Very elongated or branched shape: hard to fit into pockets and membranes.")
rule("FractionCSP3", "Fraction sp3 carbons (too low)", "below", {"*": 0.2}, "", "shape", "Almost completely flat: often poorly soluble and promiscuous.")
rule("FractionCSP3", "Fraction sp3 carbons (too high)", "above", {"*": 0.8}, "", "shape", "Very saturated and three-dimensional: can be hard to make and to bind selectively.")
rule("NumAromaticRings", "Aromatic rings", "above", {"*": 4}, "", "shape", "Many flat rings add greasiness and lower solubility.")
rule("LabuteASA", "Accessible surface area (Labute)", "above", {"*": 500}, "Å²", "surface", "A very large surface: typically high protein binding and poor absorption.")
rule("NumValenceElectrons", "Valence electrons", "above", {"*": 300}, "", "surface", "Very large electron count: usually a heavy, greasy molecule.")
# ── flexibility ──
rule("NumRotatableBonds", "Rotatable bonds", "above", {"oral": 10, "cns": 8, "iv": 15, "topical": 10, "transdermal": 8, "inhaled": 10, "inhaled_systemic": 10, "ocular": 10, "intrathecal": 15, "gut_local": 12, "buccal": 8, "rectal_vaginal": 12}, "", "flexibility",
     {"default": "Very floppy: loses energy when binding, so it binds weakly.", "oral": "Floppy molecules absorb poorly from the gut.",
      "cns": "Too floppy to pass the blood-brain barrier efficiently."})
# ── lipophilicity ──
rule("MolLogP", "Lipophilicity (logP)", "above", {"oral": 5, "cns": 4, "iv": 5, "topical": 4, "transdermal": 4, "inhaled": 3, "inhaled_systemic": 4, "ocular": 4, "intrathecal": 5, "gut_local": 5, "buccal": 4, "rectal_vaginal": 4}, "", "lipophilicity",
     {"default": "High lipophilicity can reduce aqueous solubility and increase tissue/protein binding; check formulation and exposure.",
      "topical": "Too greasy: stays trapped in the skin's outer layer instead of reaching deeper tissue.",
      "iv": "High lipophilicity can reduce aqueous/blood solubility; a formulation strategy may be needed.",
      "cns": "Too greasy: accumulates in fat and brain lipids, causing side effects."})
rule("MolLogP", "Lipophilicity (logP, too low)", "below", {"oral": 0, "cns": 1, "topical": 1, "transdermal": 1, "inhaled": 0, "inhaled_systemic": 0, "buccal": 0}, "", "lipophilicity",
     {"default": "Too water-loving: struggles to cross cell membranes.",
      "cns": "Too water-loving to cross the blood-brain barrier.", "topical": "Too water-loving to get through the skin."}, severe={"cns"})
rule("MolMR", "Molar refractivity", "above", {"*": 130}, "", "lipophilicity", "Very bulky and polarisable: usually comes with high greasiness.")
rule("BCUT2D_LOGPHI", "Highest atomic logP contribution", "above", {"*": 3.5}, "", "lipophilicity", "A very greasy patch on the surface: sticks to membranes and proteins.")
rule("FpDensityMorgan2", "Fingerprint density (Morgan r2)", "above", {"*": 2.2}, "", "lipophilicity", "Crowded, dense structure: often correlates with poor solubility.")
# ── polarity and permeability ──
rule("TPSA", "Polar surface area (TPSA)", "above", {"oral": 140, "cns": 90, "iv": 300, "topical": 150, "transdermal": 90, "inhaled": 100, "inhaled_systemic": 90, "ocular": 180, "intrathecal": 500, "gut_local": 300, "buccal": 120, "rectal_vaginal": 180}, "Å²", "polarity",
     {"default": "High polarity can limit membrane passage; formulation, transporters and route can change exposure.", "oral": "Too polar to be absorbed well from the gut.",
      "cns": "Likely to have poor passive blood-brain barrier penetration, so systemic brain action is unlikely without a delivery strategy.",
      "topical": "Too polar to penetrate the skin.", "inhaled": "Too polar to pass into lung tissue."}, severe={"cns", "transdermal", "inhaled_systemic"})
rule("TPSA", "Polar surface area (TPSA, too low for gut-restricted drug)", "below", {"gut_local": 100}, "Å²", "polarity",
     "Absorbed too easily for a drug that must stay in the intestine.")
rule("NumHDonors", "H-bond donors", "above", {"*": 10, "oral": 5, "cns": 3, "iv": 12, "topical": 5, "transdermal": 3, "inhaled": 3, "inhaled_systemic": 3, "ocular": 10, "intrathecal": 15, "gut_local": 10, "buccal": 5, "rectal_vaginal": 10}, "", "polarity",
     {"default": "Many N-H or O-H groups must shed water to enter membranes, which costs energy.", "cns": "Too many H-bond donors to cross the blood-brain barrier."}, severe={"cns"})
rule("NumHAcceptors", "H-bond acceptors", "above", {"*": 20, "oral": 10, "cns": 7, "iv": 20, "topical": 15, "transdermal": 8, "inhaled": 10, "inhaled_systemic": 10, "ocular": 15, "intrathecal": 20, "gut_local": 20, "buccal": 10, "rectal_vaginal": 15}, "", "polarity",
     "Many N or O acceptors make the molecule hold on to water: poorer membrane crossing.")
rule("NHOHCount", "N-H and O-H count", "above", {"*": 10, "oral": 5, "cns": 3, "iv": 12, "topical": 8, "transdermal": 3, "inhaled": 5, "inhaled_systemic": 3, "ocular": 10, "intrathecal": 15, "gut_local": 10, "buccal": 5, "rectal_vaginal": 10}, "", "polarity", "Many hydrogen-bonding hydroxyl and amine groups lower absorption.")
rule("NOCount", "N and O count", "above", {"*": 20, "oral": 10, "cns": 7, "iv": 20, "topical": 15, "transdermal": 10, "inhaled": 12, "inhaled_systemic": 10, "ocular": 15, "intrathecal": 20, "gut_local": 20, "buccal": 12, "rectal_vaginal": 15}, "", "polarity", "Many nitrogen and oxygen atoms: more polar, less membrane-friendly.")
rule("NumHeteroatoms", "Heteroatom count", "above", {"*": 20, "oral": 12, "cns": 9, "iv": 24, "topical": 16, "transdermal": 12, "inhaled": 15, "inhaled_systemic": 12, "ocular": 18, "intrathecal": 24, "gut_local": 24, "buccal": 15, "rectal_vaginal": 18}, "", "polarity", "Many non-carbon atoms often correlate with poorer absorption.")
rule("MaxPartialCharge", "Largest positive partial charge", "above", {"*": 0.4}, "e", "charge",
     "Strongly positive region: can trap the molecule in acidic compartments and hinder membrane crossing.")
rule("MinPartialCharge", "Largest negative partial charge", "below", {"*": -0.5}, "e", "charge",
     "Strongly negative region: likely ionised at body pH and hard to absorb.")
rule("MaxAbsPartialCharge", "Largest absolute partial charge", "above", {"*": 0.5}, "e", "charge", "Strong charge separation: likely ionised in the body.")
# ── ionisation ──
rule("fr_COO", "Carboxylic acid groups", "above", {"oral": 1, "cns": 0, "iv": 2, "topical": 1, "transdermal": 0, "inhaled": 1, "inhaled_systemic": 0, "ocular": 2, "intrathecal": 2, "gut_local": 2, "buccal": 1, "rectal_vaginal": 1}, "", "ionisation",
     {"default": "An acid group is charged at body pH.", "cns": "Acids are charged and cross the blood-brain barrier poorly."}, severe={"cns"})
rule("fr_Ar_COO", "Aromatic carboxylic acid", "above", {"oral": 0, "cns": 0, "iv": 2, "topical": 1, "transdermal": 0, "inhaled": 1, "inhaled_systemic": 0, "ocular": 2, "intrathecal": 2, "gut_local": 2, "buccal": 0, "rectal_vaginal": 1}, "", "ionisation", "Charged at body pH; tends to bind plasma proteins strongly.")
rule("fr_quatN", "Permanent positive charge (quaternary N)", "above", {"*": 0}, "", "ionisation",
     {"default": "Permanently charged: hardly crosses membranes.", "cns": "Permanently charged: cannot cross into the brain."}, severe={"oral", "cns"})
rule("fr_NH2", "Primary amine groups", "above", {"*": 2, "oral": 2, "cns": 1, "iv": 4, "topical": 2, "transdermal": 1, "inhaled": 2, "inhaled_systemic": 1, "ocular": 3, "intrathecal": 4, "gut_local": 4, "buccal": 2, "rectal_vaginal": 3}, "", "ionisation", "Several basic amines: strongly ionised, can trap in tissue and carry hERG heart-rhythm risk.")
rule("fr_phos_acid", "Phosphoric acid groups", "above", {"*": 0}, "", "ionisation", "Highly charged: very poor membrane crossing.", severe={"oral", "cns"})
rule("fr_sulfonamd", "Sulfonamide groups", "above", {"*": 1}, "", "ionisation", "Can be weakly acidic; watch for allergic reactions and altered distribution.")
# ── structural alerts (apply to every site) ──
for rd, nm, txt in [("fr_nitro", "Nitro groups", "Can be metabolised into reactive species: a common toxicity flag."),
                    ("fr_epoxide", "Epoxide rings", "Reactive ring: can alkylate DNA and proteins."),
                    ("fr_aldehyde", "Aldehyde groups", "Reactive carbonyl: can bind proteins and cause toxicity."),
                    ("fr_azide", "Azide groups", "Potentially explosive and toxic."), ("fr_isocyan", "Isocyanate groups", "Highly reactive and toxic."),
                    ("fr_isothiocyan", "Isothiocyanate groups", "Reactive: can modify proteins."), ("fr_hdrzine", "Hydrazine groups", "A known toxicity flag."),
                    ("fr_hdrzone", "Hydrazone groups", "Can hydrolyse; a toxicity flag."), ("fr_azo", "Azo groups", "Can be metabolised into toxic amines."),
                    ("fr_diazo", "Diazo groups", "Reactive and potentially toxic."), ("fr_imine", "Imine groups", "Hydrolyses in water; often unstable."),
                    ("fr_thiocyan", "Thiocyanate groups", "Potential toxicity."), ("fr_SH", "Free thiol groups", "Reacts with proteins and oxidises easily."),
                    ("fr_C_S", "Thiocarbonyl groups", "Potentially reactive sulfur group."), ("NumRadicalElectrons", "Radical electrons", "Unpaired electrons: reactive and unstable.")]:
    rule(rd, nm, "above", {"*": 0}, "", "alert", txt)
rule("fr_barbitur", "Barbiturate core", "above", {"*": 0}, "", "alert", "Sedative class with dependence risk.")
rule("fr_benzodiazepine", "Benzodiazepine core", "above", {"*": 0}, "", "alert", "Sedation and dependence risk.")
rule("fr_prisulfonamd", "Primary sulfonamide", "above", {"*": 0}, "", "alert", "Can cause sulfa-type reactions.")

# Apply a light repeated-motif threshold to every additional RDKit fr_* counter.
# These generic alerts prompt review; they do not label ordinary functional groups toxic.
_known_fragment_rules = {item[0] for item in R["list"]}
for _rd_name, _rd_fn in Descriptors.descList:
    if _rd_name.startswith("fr_") and _rd_name not in _known_fragment_rules:
        _friendly = _rd_name[3:].replace("_", " ").replace("Ar", "aromatic ")
        rule(_rd_name, "RDKit fragment: " + _friendly, "above", {"*": 3}, "count", "fragment",
             "This pattern occurs repeatedly. Review whether it changes reactivity, metabolism, selectivity or solubility in this molecule.")
        _known_fragment_rules.add(_rd_name)
# ── metabolism and clearance ──
rule("fr_para_hydroxylation", "Para-hydroxylation sites", "above", {"*": 2}, "", "metabolism", "Several spots the liver easily oxidises: may be cleared too quickly.")
rule("fr_ketone", "Ketone groups", "above", {"*": 2}, "", "metabolism", "Several ketones: possible reduction and clearance by the liver.")
rule("fr_ester", "Ester groups", "above", {"*": 2}, "", "metabolism", "Esters are broken down quickly by blood and liver enzymes.")
rule("fr_ether", "Ether groups", "above", {"*": 4}, "", "metabolism", "Many ethers add oxidation sites and shorten the drug's life.")
rule("fr_halogen", "Halogen atoms", "above", {"*": 4}, "", "metabolism", "Many halogens add greasiness and slow clearance; watch for accumulation.")
# ── synthesis complexity ──
rule("NumAtomStereoCenters", "Stereocentres", "above", {"*": 4}, "", "complexity", "Many stereocentres: hard and expensive to make as one pure form.")
rule("BertzCT", "Molecular complexity (BertzCT)", "above", {"*": 1200}, "", "complexity", "Very complex: usually hard to synthesise and scale up.")
rule("NumSpiroAtoms", "Spiro atoms", "above", {"*": 2}, "", "complexity", "Several spiro centres make synthesis harder.")
rule("NumBridgeheadAtoms", "Bridgehead atoms", "above", {"*": 2}, "", "complexity", "Several bridged ring junctions make synthesis harder.")
# ── rings ──
rule("NumAromaticHeterocycles", "Aromatic heterocycles", "above", {"*": 3}, "", "ring", "Many heteroaromatic rings: can bind unwanted targets.")
rule("NumAliphaticRings", "Aliphatic rings", "above", {"*": 4}, "", "ring", "Many saturated rings: bulky and greasy.")
rule("NumAmideBonds", "Amide bonds", "above", {"*": 4}, "", "ring", "Many amide bonds: peptide-like, often poorly absorbed and quickly broken down.")
rule("NumSaturatedHeterocycles", "Saturated heterocycles", "above", {"*": 3}, "", "ring", "Many saturated heterocycles can raise basicity and ionisation.")
# ── electronics and general drug-likeness ──
rule("MaxAbsEStateIndex", "Maximum electrotopological state", "above", {"*": 17}, "", "electronic", "A very electron-rich or electron-poor atom: a common cause of reactivity.")
rule("MinEStateIndex", "Minimum electrotopological state", "below", {"*": -4}, "", "electronic", "A strongly electron-poor atom: possible reactivity.")

RULES = R["list"]
RULE_COUNT = len(RULES)


def _values(mol: Chem.Mol) -> Dict[str, float]:
    out = {}
    for name, fn in Descriptors.descList:
        try:
            out[name] = float(fn(mol))
        except Exception:
            out[name] = float("nan")
    return out


def _limit(limits: Dict[str, float], site: str) -> Optional[float]:
    return limits[site] if site in limits else limits.get("*")


def _consequence(text, site: str) -> str:
    if isinstance(text, dict):
        return text.get(site, text.get("default", ""))
    return text


def _site_result(values: Dict[str, float], site: str) -> Dict[str, Any]:
    hits, checked = [], set()
    for rd, name, direction, limits, unit, group, cons, severe in RULES:
        lim = _limit(limits, site)
        if lim is None:
            continue
        v = values.get(rd, float("nan"))
        if v != v:
            continue
        checked.add(rd + str(lim))
        crossed = v > lim if direction == "above" else v < lim
        if not crossed:
            continue
        word = "above" if direction == "above" else "below"
        u = (" " + unit) if unit else ""
        text = f"{name} is {round(v, 2):g}{u}, {word} the limit of {lim:g}{u} for this use. {_consequence(cons, site)}"
        hits.append({"descriptor": rd, "name": name, "value": round(v, 3), "limit": lim, "direction": direction, "unit": unit,
                     "group": group, "group_label": GROUP.get(group, group), "text": text, "severe": site in severe})
    by_group: Dict[str, List[Dict[str, Any]]] = {}
    for h in hits:
        by_group.setdefault(h["group"], []).append(h)
    severe = [h["text"] for h in hits if h["severe"]]
    if not hits:
        verdict = f"No warnings for this use: no checked threshold was crossed."
    elif severe:
        verdict = ("Likely unsuitable for this intended use unless a delivery strategy changes exposure. " + " ".join(severe))
    else:
        verdict = (f"{len(hits)} warning{'s' if len(hits) != 1 else ''}. None is an absolute bar on its own, "
                   "but together they lower the chance of a usable drug for this use.")
    return {"site": site, "site_label": SITES[site], "warnings": hits, "count": len(hits), "by_group": by_group,
            "groups": [g for g in GROUP if g in by_group], "blocking": severe, "verdict": verdict,
            "descriptors_checked": len(checked), "descriptors_available": len(values)}


def warnings_for(mol: Chem.Mol, site: str = "oral") -> Dict[str, Any]:
    """Warnings for one kind of drug (site), with a verdict."""
    site = site if site in SITES else "oral"
    return _site_result(_values(Chem.RemoveHs(mol)), site)


def warnings_all(mol: Chem.Mol) -> Dict[str, Any]:
    """Warnings for every kind of drug at once, so the page can switch between them instantly."""
    values = _values(Chem.RemoveHs(mol))
    by_site = {s: _site_result(values, s) for s in SITE_ORDER}
    return {"by_site": by_site, "labels": SITES, "order": SITE_ORDER, "default": "oral",
            "summary": {s: {"count": r["count"], "suitable": not r["blocking"]} for s, r in by_site.items()},
            "descriptors_available": len(values), "rules": RULE_COUNT,
            "note": "Heuristic physicochemical screen only; not a clinical prediction or a substitute for formulation, PK or safety studies."}


def all_values(mol: Chem.Mol) -> Dict[str, Any]:
    vals = _values(Chem.RemoveHs(mol))
    return {"values": {k: (round(v, 4) if v == v else None) for k, v in vals.items()}, "count": len(vals)}
