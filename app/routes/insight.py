"""
routes/insight.py — Chemistry Guide: pocket chemistry, molecule profile, match-maker, interaction map.

Everything here is rule-of-thumb chemistry meant for learning — it is NOT a validated predictor.
"""
import asyncio
import logging
from pathlib import Path
from typing import List

import numpy as np
from fastapi import APIRouter, HTTPException
from pydantic import BaseModel
from rdkit import Chem
from rdkit.Chem import Crippen, Descriptors, Lipinski, rdMolDescriptors

from .. import docking as dock_engine
from ..config import RECEPTORS_DIR
from ..lipinski import lipinski_filter
from ..synth import synth_summary
from ..receptor_prep import RECEPTORS

logger = logging.getLogger("medchem.insight")
router = APIRouter()

# ─────────────────────────── residue chemistry ───────────────────────────
RES_CLASS = {r: "hydrophobic" for r in "ALA VAL LEU ILE MET PRO".split()}
RES_CLASS.update({r: "aromatic" for r in "PHE TRP TYR".split()})
RES_CLASS.update({r: "polar" for r in "SER THR ASN GLN CYS HIS".split()})
RES_CLASS.update({"ASP": "acidic", "GLU": "acidic", "LYS": "basic", "ARG": "basic", "GLY": "small"})
BACKBONE = {"N", "CA", "C", "O", "OXT"}
METALS = {"FE", "ZN", "MG", "MN", "CU", "NI", "CO", "CA"}
COFACTOR_NOTE = {
    "FAD": "FAD — the yellow flavin cofactor that does the chemistry (it oxidises amines). Flat aromatic rings and "
           "propargyl (–C≡C) groups next to it are classic MAO-inhibitor features.",
    "HEM": "Heme — the iron-containing reaction centre. Nitrogen atoms with a free lone pair (imidazole, pyridine) "
           "can bind the iron directly — that is how azole drugs block CYP enzymes.",
    "NAP": "NADP — a nicotinamide cofactor (carries hydride)", "NDP": "NADPH — a nicotinamide cofactor",
    "NAD": "NAD — a nicotinamide cofactor",
}
LANDMARKS = {   # target PDB id → (residue, number, note); shown only if that exact residue lines the box
    "2Z5Y": [("TYR", 407, "Aromatic-cage tyrosine — flat rings and amines sit against it"),
             ("TYR", 444, "Aromatic-cage tyrosine"), ("PHE", 208, "Greasy wall of the substrate cavity"),
             ("ILE", 335, "Greasy residue that helps pick small, flat ligands")],
    "2V5Z": [("TYR", 398, "Aromatic-cage tyrosine — flat rings and amines sit against it"),
             ("TYR", 435, "Aromatic-cage tyrosine"),
             ("ILE", 199, "'Gate' residue between the entrance and substrate cavities")],
    "4JIT": [("TYR", 71, "'Flap' tyrosine — the lid that closes over the ligand")],
    "1W0G": [("PHE", 304, "I-helix phenylalanine near the heme — part of the greasy 'Phe cluster'")],
}

# friendly chemistry features (also used by the Drug-Twin similarity card)
FEATS = [  # key, emoji, label, SMARTS, tip
    ("arom", "🔵", "Aromatic ring", "a", "Flat ring with a cloud of π electrons — stacks on Phe/Tyr/Trp rings."),
    ("acid", "🔴", "Carboxylic acid (−)", "[CX3](=O)[OX2H1,OX1-]", "Negative at body pH — looks for Lys/Arg (+)."),
    ("base", "🔷", "Basic amine (+)",
     "[NX3;!a;!$(N-C=[O,S,N]);!$(N-a);!$(N-S);!$(N-[#7,#8]);!$(N-C#N)]", "Positive at body pH — looks for Asp/Glu (−)."),
    ("amidine", "🔷", "Amidine / guanidine (+)", "[NX3;!a][CX3;!a]=[NX2;!a]", "Strongly positive — loves Asp/Glu."),
    ("oh", "💧", "Hydroxyl (–OH)", "[OX2H][CX4,c]", "H-bond donor AND acceptor."),
    ("amide", "🧷", "Amide", "C(=O)[NX3]", "N–H donates, C=O accepts H-bonds; rigid and flat."),
    ("carbonyl", "⚔️", "Carbonyl (C=O)", "[CX3]=[OX1]", "H-bond acceptor."),
    ("ether", "🌉", "Ether (–O–)", "[OD2]([#6])[#6]", "H-bond acceptor."),
    ("halogen", "🟢", "Halogen (F/Cl/Br/I)", "[F,Cl,Br,I]", "Greasy; fits tight oily pockets (Cl/Br can form halogen bonds)."),
    ("sulfur", "🟡", "Sulfur", "[#16]", "Sulfonamides/thioethers — polar or greasy depending on type."),
    ("azole", "🎯", "Ring N with free lone pair", "[nX2;r5,r6]", "Can grab a metal — e.g. the iron of heme in CYP enzymes."),
    ("alkyne", "🪝", "Alkyne (C≡C)", "[CX2]#[CX2]", "Propargyl groups can bond covalently to FAD (selegiline-style)."),
    ("nitrile", "🧪", "Nitrile (C≡N)", "[CX2]#[NX1]", "Small, polar H-bond acceptor."),
    ("ring", "⭕", "Non-aromatic ring", "[R;!a]", "Gives 3-D shape and rigidity."),
]
_FEAT_PATT = [(k, e, l, Chem.MolFromSmarts(s), t) for k, e, l, s, t in FEATS]
FEAT_LABEL = {k: (e, l, t) for k, e, l, _, t in _FEAT_PATT}


def feature_keys(mol) -> set:
    return {k for k, _, _, p, _ in _FEAT_PATT if p is not None and mol.HasSubstructMatch(p)}


# ─────────────────────────── PDB helpers ───────────────────────────
def _parse_pdb_lines(lines) -> list:
    atoms = []
    for ln in lines:
        if not ln.startswith(("ATOM", "HETATM")) or ln[16:17] not in (" ", "A", "1"):
            continue
        res = ln[17:20].strip()
        if res in ("HOH", "WAT", "DOD"):
            continue
        try:
            xyz = (float(ln[30:38]), float(ln[38:46]), float(ln[46:54]))
        except ValueError:
            continue
        name = ln[12:16].strip()
        el = (ln[76:78].strip() or name.lstrip("0123456789")[:1]).upper()
        if el == "H":
            continue
        try:
            resi = int(ln[22:26])
        except ValueError:
            resi = 0
        atoms.append({"name": name, "res": res, "chain": ln[21:22].strip() or "A", "resi": resi,
                      "icode": ln[26:27].strip(), "xyz": xyz, "el": el, "het": ln.startswith("HETATM")})
    return atoms


_CACHE = {}


def _receptor_atoms(tid: str):
    info = RECEPTORS.get(tid)
    if not info:
        raise HTTPException(status_code=404, detail=f"Unknown target: {tid}")
    p = Path(info["pdb_path"]) if "pdb_path" in info else RECEPTORS_DIR / f"{tid}_{info['pdb_id']}_clean.pdb"
    if not p.exists():
        raise HTTPException(status_code=404, detail=f"Receptor PDB not found for {tid}")
    key = (str(p), p.stat().st_mtime)
    if key not in _CACHE:
        _CACHE.clear()
        _CACHE[key] = _parse_pdb_lines(p.read_text(encoding="utf-8", errors="ignore").splitlines())
    return _CACHE[key], str(info.get("pdb_id", "")).upper()


def _box(tid: str):
    b = dock_engine.BOX_PARAMS.get(tid)
    if not b:
        raise HTTPException(status_code=404, detail=f"No grid box for {tid}")
    return b["center"], b["size"]


def _group(atoms):
    groups = {}
    for a in atoms:
        groups.setdefault((a["chain"], a["resi"], a["icode"], a["res"]), []).append(a)
    return groups


# ─────────────────────────── 1) pocket chemistry ───────────────────────────
def analyse_pocket(tid: str) -> dict:
    atoms, pdb_id = _receptor_atoms(tid)
    center, size = _box(tid)
    margin = 2.0
    half = [size[i] / 2 + margin for i in range(3)]
    groups = _group(atoms)

    def inside(a):
        return all(abs(a["xyz"][i] - center[i]) <= half[i] for i in range(3))

    residues, cofactors = [], []
    for (chain, resi, icode, res), ats in groups.items():
        cls = RES_CLASS.get(res)
        if cls is None:                                   # cofactor / ion / other hetero group
            if not any(inside(a) for a in ats):
                continue
            if len(ats) == 1 and ats[0]["el"] in METALS:
                cls = "metal"
            elif len(ats) >= 6:
                cls = "cofactor"
            else:
                continue
            if cls == "cofactor" and res not in cofactors:
                cofactors.append(res)
        else:
            side = [a for a in ats if a["name"] not in BACKBONE] or ats
            if not any(inside(a) for a in (side if res != "GLY" else ats)):
                continue
        ca = next((a for a in ats if a["name"] == "CA"), ats[0])
        side = [a for a in ats if a["name"] not in BACKBONE] or ats
        cen = np.mean([a["xyz"] for a in side], axis=0)
        residues.append({"chain": chain, "resi": resi, "res": res, "cls": cls, "xyz": [round(v, 2) for v in ca["xyz"]],
                         "dist": round(float(np.linalg.norm(cen - np.array(center))), 1), "note": ""})
    residues.sort(key=lambda r: r["dist"])

    # landmarks: known residues for built-in targets + the aspartyl-protease "DTG/DSG" signature
    notes = {(n[0], n[1]): n[2] for n in LANDMARKS.get(pdb_id, [])}
    by_chain = {}
    for (chain, resi, icode, res) in groups:
        by_chain.setdefault(chain, {})[resi] = res
    for r in residues:
        if (r["res"], r["resi"]) in notes:
            r["note"] = notes[(r["res"], r["resi"])]
        elif r["res"] == "ASP":
            nxt = by_chain.get(r["chain"], {})
            if nxt.get(r["resi"] + 1) in ("THR", "SER") and nxt.get(r["resi"] + 2) == "GLY":
                r["note"] = ("Catalytic Asp (D-T/S-G motif) — one half of the 'acid–acid scissors' that cut the substrate. "
                             "Positive groups (protonated amines, amidines) bind here.")
        if r["cls"] == "cofactor" and r["res"] in COFACTOR_NOTE:
            r["note"] = COFACTOR_NOTE[r["res"]]

    std = [r for r in residues if r["cls"] in ("hydrophobic", "aromatic", "polar", "acidic", "basic", "small")]
    counts = {c: sum(1 for r in residues if r["cls"] == c)
              for c in ("hydrophobic", "aromatic", "polar", "acidic", "basic", "small", "cofactor", "metal")}
    n = max(1, len(std))
    greasy = round(100 * (counts["hydrophobic"] + counts["aromatic"]) / n)
    net = counts["basic"] - counts["acidic"]
    if greasy >= 55:
        character = "mostly greasy (oily) — greasy molecules feel at home here"
    elif greasy <= 30:
        character = "mostly polar/charged — it likes water-friendly molecules with H-bonding groups"
    else:
        character = "a mix of greasy and polar zones"

    hints = []
    if counts["acidic"]:
        hints.append(f"{counts['acidic']} acidic residue(s) (−): a basic amine or amidine (+) on your molecule can form a salt bridge.")
    if counts["basic"]:
        hints.append(f"{counts['basic']} basic residue(s) (+): an acid group (−) on your molecule would be attracted.")
    if counts["aromatic"]:
        hints.append(f"{counts['aromatic']} aromatic residue(s): flat aromatic rings on your molecule can stack against them.")
    if counts["polar"]:
        hints.append(f"{counts['polar']} polar residue(s): –OH, –NH and C=O groups on your molecule can H-bond with them.")
    if counts["hydrophobic"] >= 4:
        hints.append("Lots of greasy side chains: carbon-rich, low-polarity parts of your molecule will be happy here.")
    if "HEM" in cofactors:
        hints.append("Heme in the box: a ring N with a free lone pair (imidazole, pyridine) can coordinate the iron.")
    if "FAD" in cofactors:
        hints.append("FAD in the box: flat aromatic rings and propargyl (–C≡C) groups are classic partners.")
    return {"target_id": tid, "box": {"center": list(center), "size": list(size)}, "n_residues": len(residues),
            "counts": counts, "greasy_pct": greasy, "net_charge": net, "character": character,
            "residues": residues, "cofactors": cofactors, "hints": hints}


# ─────────────────────────── 2) molecule profile ───────────────────────────
def profile_smiles(smiles: str) -> dict:
    mol = Chem.MolFromSmiles((smiles or "").strip())
    if mol is None:
        raise ValueError("That SMILES is not valid")
    keys = feature_keys(mol)
    logp = Crippen.MolLogP(mol)
    ha = mol.GetNumHeavyAtoms()
    return {
        "smiles": smiles.strip(), "heavy_atoms": ha, "mw": round(Descriptors.MolWt(mol), 1),
        "logp": round(logp, 2), "tpsa": round(rdMolDescriptors.CalcTPSA(mol), 1),
        "hbd": Lipinski.NumHDonors(mol), "hba": Lipinski.NumHAcceptors(mol),
        "rot_bonds": Lipinski.NumRotatableBonds(mol), "arom_rings": rdMolDescriptors.CalcNumAromaticRings(mol),
        "greasiness": "water-loving" if logp < 1 else ("balanced" if logp <= 3 else "greasy"),
        "size_label": "tiny fragment" if ha < 12 else ("small" if ha < 20 else ("drug-sized" if ha <= 38 else "large")),
        "lipinski": lipinski_filter(mol, "standard"),
        "synth": synth_summary(mol),
        "feature_keys": sorted(keys),
        "features": [{"key": k, "emoji": FEAT_LABEL[k][0], "label": FEAT_LABEL[k][1], "tip": FEAT_LABEL[k][2]}
                     for k, *_ in FEATS if k in keys],
    }


# ─────────────────────────── 3) match-maker ───────────────────────────
def match_score(pocket: dict, prof: dict) -> dict:
    c, f = pocket["counts"], set(prof["feature_keys"])
    pos, neg = bool(f & {"base", "amidine"}), "acid" in f
    pts, reasons = 50, []

    def add(icon, text, d):
        nonlocal pts
        pts += d
        reasons.append({"icon": icon, "text": text})

    if c["acidic"] and pos:
        add("✅", f"Opposites attract: your molecule has a positive group (+) and the pocket has {c['acidic']} acidic residue(s) (−).", 14)
    elif c["acidic"] and neg and not c["basic"]:
        add("⚠️", "Like charges repel: your molecule is acidic (−) and the pocket is rich in acidic residues (−).", -12)
    if c["basic"] and neg:
        add("✅", f"Opposites attract: your acid group (−) can pair with the pocket's {c['basic']} basic residue(s) (+).", 12)
    elif c["basic"] and pos and not c["acidic"]:
        add("⚠️", "Like charges repel: your molecule is positive (+) and the pocket has basic residues (+) but no acidic partner.", -10)
    if c["aromatic"] and prof["arom_rings"] >= 1:
        add("✅", "Your aromatic ring can stack against the pocket's aromatic residues (π-stacking).", 8)
    elif c["aromatic"] >= 2 and prof["arom_rings"] == 0:
        add("⚠️", "The pocket has several aromatic residues but your molecule has no aromatic ring to stack on them.", -6)
    n_pol = c["polar"] + c["acidic"] + c["basic"]
    if n_pol >= 3 and prof["hbd"] + prof["hba"] >= 3:
        add("✅", "The pocket has polar spots and your molecule has H-bond donors/acceptors to match.", 8)
    elif n_pol >= 4 and prof["hbd"] + prof["hba"] <= 1:
        add("⚠️", "The pocket is full of polar/charged residues but your molecule has almost nothing to H-bond with.", -6)
    if pocket["greasy_pct"] >= 55:
        if prof["logp"] >= 2:
            add("✅", "Greasy pocket, greasy molecule (logP ≥ 2): oil likes oil.", 8)
        elif prof["logp"] < 1:
            add("⚠️", "Greasy pocket but a water-loving molecule (logP < 1): it would rather stay in the water.", -10)
    elif pocket["greasy_pct"] <= 30 and prof["logp"] > 4.5:
        add("⚠️", "Polar pocket but a very greasy molecule (logP > 4.5).", -6)
    ha = prof["heavy_atoms"]
    if ha < 12:
        add("⚠️", f"Tiny molecule ({ha} heavy atoms): few contacts, so usually weak binding.", -14)
    elif ha < 20:
        add("💡", f"Smallish molecule ({ha} heavy atoms): bigger molecules make more contacts and usually score better.", -4)
    elif ha <= 38:
        add("✅", f"Good size ({ha} heavy atoms): large enough for many contacts, small enough to fit.", 10)
    elif ha > 45:
        add("⚠️", f"Very large molecule ({ha} heavy atoms): may not fit the box and is very floppy.", -8)
    lip = prof.get("lipinski")
    if lip and lip["violations"] > 1:
        add("⚠️", f"Misses {lip['violations']} Lipinski Rule-of-Five limits ({'; '.join(lip['violations_detail'])}) — hard to turn into an oral pill.", -8)
    elif lip and lip["violations"] == 1:
        add("💡", f"Misses one Rule-of-Five limit ({lip['violations_detail'][0]}) — still acceptable, but keep an eye on it.", -2)
    elif lip:
        add("✅", "Passes Lipinski's Rule of Five — the right size and polarity for an oral pill.", 3)
    syn = prof.get("synth")
    if syn and syn.get("available"):
        if syn["tier"] == "very_hard":
            add("⚠️", f"Very hard to synthesize (SA score {syn['sa']}): a great pose is worth little if nobody can make the molecule.", -10)
        elif syn["tier"] == "hard":
            add("💡", f"Challenging to make (SA score {syn['sa']}): simplify stereocentres, rings or size if you can.", -4)
        elif syn["tier"] == "easy":
            add("✅", f"Easy to synthesize (SA score {syn['sa']}): a chemist could make this quickly.", 3)
    if prof["rot_bonds"] > 10:
        add("⚠️", f"{prof['rot_bonds']} rotatable bonds: a floppy molecule pays an entropy penalty to freeze into the pocket.", -5)
    if "HEM" in pocket["cofactors"]:
        if "azole" in f:
            add("✅", "Heme iron in the box and your molecule has a ring N with a free lone pair — it can bind the iron.", 14)
        else:
            add("💡", "Tip: heme is in the box — a pyridine/imidazole-type N could bind its iron.", 0)
    if "FAD" in pocket["cofactors"]:
        if "alkyne" in f:
            add("✅", "FAD in the box and your molecule has a propargyl alkyne — the selegiline-style trick.", 12)
        elif prof["arom_rings"]:
            add("💡", "FAD in the box: your flat aromatic ring can sit against the flavin.", 4)
    score = max(5, min(95, pts))
    return {"score": score, "category": "strong" if score >= 68 else ("medium" if score >= 45 else "weak"),
            "verdict": "Looks promising" if score >= 68 else ("Could go either way" if score >= 45 else "Looks like a poor fit"),
            "reasons": reasons}


# ─────────────────────────── 4) interaction map ───────────────────────────
CHARGED = {("ASP", "OD1"): "−", ("ASP", "OD2"): "−", ("GLU", "OE1"): "−", ("GLU", "OE2"): "−",
           ("LYS", "NZ"): "+", ("ARG", "NE"): "+", ("ARG", "NH1"): "+", ("ARG", "NH2"): "+"}
HYD_RES = {"ALA", "VAL", "LEU", "ILE", "MET", "PRO", "PHE", "TRP", "TYR"}
RING_ATOMS = {"PHE": [["CG", "CD1", "CD2", "CE1", "CE2", "CZ"]], "TYR": [["CG", "CD1", "CD2", "CE1", "CE2", "CZ"]],
              "HIS": [["CG", "ND1", "CD2", "CE1", "NE2"]],
              "TRP": [["CG", "CD1", "NE1", "CE2", "CD2"], ["CD2", "CE2", "CE3", "CZ2", "CZ3", "CH2"]]}


def _plane(pts):
    c = pts.mean(axis=0)
    u, s, vt = np.linalg.svd(pts - c)
    n = vt[2]
    return c, n, float(np.abs((pts - c) @ n).max())


def analyse_interactions(tid: str, docked_pdb: str) -> dict:
    prot, _ = _receptor_atoms(tid)
    lig = _parse_pdb_lines(docked_pdb.splitlines())
    if not lig:
        raise ValueError("The docked pose has no atoms")
    P = np.array([a["xyz"] for a in prot])
    L = np.array([a["xyz"] for a in lig])
    D = np.linalg.norm(L[:, None, :] - P[None, :, :], axis=2)
    items, seen = [], set()

    def push(kind, a, detail, lxyz, pxyz, dist):
        items.append({"type": kind, "label": f"{a['res'].title()}{a['resi']}", "res": a["res"], "chain": a["chain"],
                      "resi": a["resi"], "detail": detail, "lig": [round(float(v), 2) for v in lxyz],
                      "prot": [round(float(v), 2) for v in pxyz], "dist": round(float(dist), 2)})

    # H-bonds / charge pairs / metal contacts
    for i, la in enumerate(lig):
        for j in np.where(D[i] <= 4.0)[0]:
            pa, d = prot[j], D[i][j]
            ch = CHARGED.get((pa["res"], pa["name"]))
            if pa["el"] in METALS and pa["het"] and d <= 2.8 and la["el"] in ("N", "O", "S"):
                k = ("metal", pa["resi"], i)
                if k not in seen:
                    seen.add(k); push("metal", pa, f"{la['el']} binds the metal at {d:.1f} Å", L[i], P[j], d)
            elif la["el"] in ("N", "O") and pa["el"] in ("N", "O") and d <= 4.0:
                paired = ch and ((ch == "−" and la["el"] == "N") or (ch == "+" and la["el"] == "O"))
                if paired and d <= 4.0:
                    k = ("ionic", pa["chain"], pa["resi"])
                    if k not in seen:
                        seen.add(k); push("ionic", pa, f"charge pair ({ch} residue), {d:.1f} Å", L[i], P[j], d)
                if d <= 3.5 and d >= 2.4:
                    k = ("hbond", pa["chain"], pa["resi"], i)
                    if k not in seen:
                        seen.add(k); push("hbond", pa, f"{la['el']}···{pa['name']} {d:.1f} Å", L[i], P[j], d)

    # greasy contacts (≥3 carbon–carbon pairs within 4 Å)
    per = {}
    for i, la in enumerate(lig):
        if la["el"] != "C":
            continue
        for j in np.where(D[i] <= 4.0)[0]:
            pa = prot[j]
            if pa["el"] == "C" and pa["res"] in HYD_RES and pa["name"] not in BACKBONE:
                per.setdefault((pa["chain"], pa["resi"], pa["res"]), []).append((i, j))
    for (chain, resi, res), pairs in per.items():
        if len(pairs) >= 3:
            li = np.mean([L[i] for i, _ in pairs], axis=0)
            pj = np.mean([P[j] for _, j in pairs], axis=0)
            push("hydrophobic", {"res": res, "resi": resi, "chain": chain}, f"{len(pairs)} greasy contacts", li, pj,
                 np.linalg.norm(li - pj))

    # ring stacking: planar ligand rings vs aromatic residue rings
    try:
        m = Chem.MolFromPDBBlock("\n".join(l for l in docked_pdb.splitlines() if l.startswith(("ATOM", "HETATM"))),
                                 sanitize=False, removeHs=True)
        m.UpdatePropertyCache(strict=False)
        Chem.FastFindRings(m)
        conf = m.GetConformer()
        lig_rings = []
        for ring in m.GetRingInfo().AtomRings():
            if len(ring) in (5, 6) and all(m.GetAtomWithIdx(k).GetSymbol() in ("C", "N") for k in ring):
                pts = np.array([list(conf.GetAtomPosition(k)) for k in ring])
                c, nrm, dev = _plane(pts)
                if dev < 0.12:
                    lig_rings.append((c, nrm))
        by_res = _group(prot)
        for (chain, resi, icode, res), ats in by_res.items():
            for names in RING_ATOMS.get(res, []):
                pts = np.array([a["xyz"] for a in ats if a["name"] in names])
                if len(pts) < 5:
                    continue
                pc, pn, _ = _plane(pts)
                for lc, ln_ in lig_rings:
                    d = float(np.linalg.norm(lc - pc))
                    ang = float(np.degrees(np.arccos(min(1.0, abs(float(ln_ @ pn))))))
                    if d <= 5.5 and (ang <= 35 or ang >= 55):
                        k = ("pi", chain, resi)
                        if k not in seen:
                            seen.add(k)
                            push("pi", {"res": res, "resi": resi, "chain": chain},
                                 ("face-to-face stack" if ang <= 35 else "T-shaped stack") + f", {d:.1f} Å", lc, pc, d)
    except Exception as exc:
        logger.info("ring-stacking check skipped: %s", exc)

    near = (D <= 4.5).sum(axis=1)
    buried = round(100 * float((near >= 8).mean()))
    kinds = {k: sum(1 for x in items if x["type"] == k) for k in ("hbond", "ionic", "hydrophobic", "pi", "metal")}
    items.sort(key=lambda x: ("metal", "ionic", "hbond", "pi", "hydrophobic").index(x["type"]))
    return {"items": items, "counts": kinds, "buried_pct": buried, "n_ligand_atoms": len(lig),
            "n_contact_residues": len({(x["chain"], x["resi"]) for x in items})}


# ─────────────────────────── API ───────────────────────────
class SmilesReq(BaseModel):
    smiles: str


class MatchReq(BaseModel):
    target_id: str
    smiles: str


class InterReq(BaseModel):
    target_id: str
    docked_pdb: str


def _bad(exc):
    raise HTTPException(status_code=422, detail=str(exc))


@router.get("/api/insight/pocket/{target_id}")
async def api_pocket(target_id: str):
    return await asyncio.to_thread(analyse_pocket, target_id)


@router.post("/api/insight/ligand")
async def api_ligand(req: SmilesReq):
    try:
        return await asyncio.to_thread(profile_smiles, req.smiles)
    except ValueError as exc:
        _bad(exc)


@router.post("/api/insight/match")
async def api_match(req: MatchReq):
    try:
        def work():
            return match_score(analyse_pocket(req.target_id), profile_smiles(req.smiles))
        return await asyncio.to_thread(work)
    except ValueError as exc:
        _bad(exc)


@router.post("/api/insight/interactions")
async def api_interactions(req: InterReq):
    try:
        return await asyncio.to_thread(analyse_interactions, req.target_id, req.docked_pdb)
    except ValueError as exc:
        _bad(exc)


# ─────────────────────────── 5) ligand preview + protonated receptor (viewer) ───────────────────────────
from fastapi.responses import FileResponse   # noqa: E402


class PreviewReq(BaseModel):
    target_id: str
    smiles: str


@router.post("/api/insight/preview")
async def api_preview(req: PreviewReq):
    """3-D conformer of the molecule, centred in the grid box, so it can be viewed BEFORE docking."""
    def work():
        from ..chem import smiles_to_mol3d
        mol = smiles_to_mol3d(req.smiles.strip())
        center, size = _box(req.target_id)
        conf = mol.GetConformer()
        xyz = conf.GetPositions()
        shift = np.array(center) - xyz.mean(axis=0)
        for i in range(mol.GetNumAtoms()):
            conf.SetAtomPosition(i, (xyz[i] + shift).tolist())
        heavy = np.array([xyz[i] + shift for i, a in enumerate(mol.GetAtoms()) if a.GetAtomicNum() > 1])
        length = float(np.linalg.norm(heavy[:, None] - heavy[None], axis=2).max()) if len(heavy) > 1 else 0.0
        return {"sdf": Chem.MolToMolBlock(mol), "heavy_atoms": int(len(heavy)), "length": round(length, 1),
                "box_min_edge": round(float(min(size)), 1), "fits": length <= float(min(size))}
    try:
        return await asyncio.to_thread(work)
    except HTTPException:
        raise
    except Exception as exc:
        raise HTTPException(status_code=422, detail=str(exc) or "Could not build a 3-D model of this molecule")


@router.get("/api/receptor/{target_id}/pdb_h")
async def receptor_pdb_h(target_id: str):
    """The protonated (hydrogens added) receptor file, for the 'Hydrogens' display option."""
    info = RECEPTORS.get(target_id)
    if not info:
        raise HTTPException(status_code=404, detail=f"Unknown target: {target_id}")
    clean = Path(info["pdb_path"]) if "pdb_path" in info else RECEPTORS_DIR / f"{target_id}_{info['pdb_id']}_clean.pdb"
    hp = clean.with_name(clean.name.replace("_clean.pdb", "_H.pdb"))
    if not hp.exists():
        raise HTTPException(status_code=404, detail="Protonated receptor file not found")
    return FileResponse(str(hp), media_type="text/plain")
