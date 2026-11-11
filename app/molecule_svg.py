"""Interactive RDKit SVG rendering for molecular similarity results."""
import html
import json
import re
from collections import Counter

from rdkit import Chem


_FUNCTIONAL_GROUPS = (
    ("Carboxylic acid (–COOH)", "[CX3](=O)[OX2H1]"),
    ("Carboxylate (–COO⁻)", "[CX3](=O)[O-]"),
    ("Carbamate", "[NX3][CX3](=O)[OX2][#6]"),
    ("Urea", "[NX3][CX3](=O)[NX3]"),
    ("Amide", "[CX3](=O)[NX3]"),
    ("Ester", "[CX3](=O)[OX2H0][#6]"),
    ("Aldehyde (–CHO)", "[CX3H1](=O)[#6]"),
    ("Ketone (C=O)", "[#6][CX3](=O)[#6]"),
    ("Phenol (Ar–OH)", "[c][OX2H1]"),
    ("Alcohol (–OH)", "[#6][OX2H]"),
    ("Ether (C–O–C)", "[#6][OD2][#6]"),
    ("Sulfonamide", "[SX4](=O)(=O)[NX3]"),
    ("Sulfonic acid", "[SX4](=O)(=O)[OX2H1]"),
    ("Sulfone", "[#6][SX4](=O)(=O)[#6]"),
    ("Sulfoxide", "[#6][SX3](=O)[#6]"),
    ("Thioether (C–S–C)", "[#6][SX2;!a][#6]"),
    ("Nitro group (–NO₂)", "[N+](=O)[O-]"),
    ("Nitrile (–C≡N)", "[#6]#[NX1]"),
    ("Amidine / guanidine", "[NX3;!a][CX3;!a]=[NX2;!a]"),
    ("Primary amine", "[NX3;H2;!$(N-C=O)]"),
    ("Secondary amine", "[NX3;H1;!$(N-C=O)]"),
    ("Tertiary amine", "[NX3;H0;!$(N-C=O)]"),
    ("Thiol (–SH)", "[SX2H]"),
    ("Azole ring nitrogen", "[nX2;r5,r6]"),
    ("Halogen (F/Cl/Br/I)", "[F,Cl,Br,I]"),
)


def _atom_interaction_info(mol):
    """Map atoms to a concise functional-group label for the interactive SVG."""
    info = [None] * mol.GetNumAtoms()
    for label, smarts in _FUNCTIONAL_GROUPS:
        pattern = Chem.MolFromSmarts(smarts)
        if pattern is None:
            continue
        for match in mol.GetSubstructMatches(pattern):
            group_atoms = list(match)
            try:
                group_smiles = Chem.MolFragmentToSmiles(
                    mol, atomsToUse=group_atoms, canonical=True, isomericSmiles=True)
            except Exception:
                group_smiles = ""
            for index in group_atoms:
                if info[index] is None:
                    atom = mol.GetAtomWithIdx(index)
                    info[index] = {
                        "index": index,
                        "element": atom.GetSymbol(),
                        "label": label,
                        "group_smiles": group_smiles,
                        "group_atom_indices": group_atoms,
                    }
    ring_info = {}
    for ring in mol.GetRingInfo().AtomRings():
        atoms = sorted(ring)
        if not atoms:
            continue
        label = "Aromatic ring" if all(mol.GetAtomWithIdx(i).GetIsAromatic() for i in atoms) else "Non-aromatic ring"
        try:
            ring_smiles = Chem.MolFragmentToSmiles(mol, atomsToUse=atoms, canonical=True, isomericSmiles=True)
        except Exception:
            ring_smiles = ""
        for index in atoms:
            if index not in ring_info or len(atoms) > len(ring_info[index][1]):
                ring_info[index] = (label, atoms, ring_smiles)

    for index, item in enumerate(info):
        if item is not None:
            continue
        atom = mol.GetAtomWithIdx(index)
        if index in ring_info:
            label, group_atoms, group_smiles = ring_info[index]
        elif atom.GetIsAromatic():
            label = f"Aromatic {atom.GetSymbol()}"
            group_atoms, group_smiles = [index], ""
        elif atom.IsInRing():
            label = f"{atom.GetSymbol()} in a ring"
            group_atoms, group_smiles = [index], ""
        else:
            label = f"{atom.GetSymbol()} atom"
            group_atoms, group_smiles = [index], ""
        info[index] = {
            "index": index,
            "element": atom.GetSymbol(),
            "label": label,
            "group_smiles": group_smiles,
            "group_atom_indices": group_atoms,
        }
    return info


def _interactive_svg(svg, mol):
    atom_info = _atom_interaction_info(mol)
    graphic_tag = re.compile(r"<(?:path|text|line|polyline|polygon|circle|ellipse)\b[^>]*>")

    def add_atom_data(match):
        tag = match.group(0)
        class_match = re.search(r"\bclass=(['\"])(.*?)\1", tag)
        if not class_match:
            return tag
        indices = []
        for raw_index in re.findall(r"(?:^|\s)atom-(\d+)(?=\s|$)", class_match.group(2)):
            index = int(raw_index)
            if index < len(atom_info) and index not in indices:
                indices.append(index)
        if not indices:
            return tag
        metadata = html.escape(json.dumps([atom_info[i] for i in indices], ensure_ascii=False), quote=True)
        attributes = f' data-mol-atoms="{metadata}" tabindex="0" role="button" focusable="true"'
        return tag[:-2] + attributes + "/>" if tag.endswith("/>") else tag[:-1] + attributes + ">"

    return graphic_tag.sub(add_atom_data, svg)


def _molecule_groups(mol):
    """Return named functional groups and rings as atom-index sets."""
    functional = []
    for priority, (label, smarts) in enumerate(_FUNCTIONAL_GROUPS):
        pattern = Chem.MolFromSmarts(smarts)
        if pattern is None:
            continue
        for match in mol.GetSubstructMatches(pattern, uniquify=True):
            if len(match) > 1 or label in ("Azole ring nitrogen", "Halogen (F/Cl/Br/I)"):
                functional.append((priority, label, frozenset(match)))

    # Keep the most specific match when generic SMARTS overlap a named group
    # (for example, an ester also matches the broad C-O-C ether pattern).
    groups = set()
    for priority, label, atoms in functional:
        shadowed = any(
            other_priority < priority
            and atoms.issubset(other_atoms)
            and (atoms != other_atoms or other_label != label)
            for other_priority, other_label, other_atoms in functional
        )
        if not shadowed:
            groups.add((label, atoms))
    for ring in mol.GetRingInfo().AtomRings():
        atoms = frozenset(ring)
        if atoms:
            ring_type = "Aromatic" if all(mol.GetAtomWithIdx(i).GetIsAromatic() for i in atoms) else "Non-aromatic"
            label = f"{ring_type} ring ({len(atoms)} atoms)"
            groups.add((label, atoms))
    return groups


def compare_molecular_features(query, candidate):
    """Summarize shared and side-specific functional groups and ring types."""
    query_counts = Counter(label for label, _ in _molecule_groups(query))
    candidate_counts = Counter(label for label, _ in _molecule_groups(candidate))
    labels = sorted(set(query_counts) | set(candidate_counts))

    shared = [{"label": label, "count": min(query_counts[label], candidate_counts[label])}
              for label in labels if query_counts[label] and candidate_counts[label]]
    only_query = [{"label": label, "count": query_counts[label] - candidate_counts[label]}
                  for label in labels if query_counts[label] > candidate_counts[label]]
    only_candidate = [{"label": label, "count": candidate_counts[label] - query_counts[label]}
                      for label in labels if candidate_counts[label] > query_counts[label]]
    return {"shared": shared, "only_query": only_query, "only_candidate": only_candidate}


def _molecule_groups(mol):
    """Return named functional groups and rings as atom-index sets."""
    functional = []
    for priority, (label, smarts) in enumerate(_FUNCTIONAL_GROUPS):
        pattern = Chem.MolFromSmarts(smarts)
        if pattern is None:
            continue
        for match in mol.GetSubstructMatches(pattern, uniquify=True):
            if len(match) > 1 or label in ("Azole ring nitrogen", "Halogen (F/Cl/Br/I)"):
                functional.append((priority, label, frozenset(match)))
    groups = set()
    for priority, label, atoms in functional:
        shadowed = any(
            other_priority < priority
            and atoms.issubset(other_atoms)
            and (atoms != other_atoms or other_label != label)
            for other_priority, other_label, other_atoms in functional
        )
        if not shadowed:
            groups.add((label, atoms))
    for ring in mol.GetRingInfo().AtomRings():
        atoms = frozenset(ring)
        if atoms:
            ring_type = "Aromatic" if all(mol.GetAtomWithIdx(i).GetIsAromatic() for i in atoms) else "Non-aromatic"
            label = f"{ring_type} ring ({len(atoms)} atoms)"
            groups.add((label, atoms))
    return groups


def compare_molecular_features(query, candidate):
    """Summarize shared and side-specific functional groups and ring types."""
    query_counts = Counter(label for label, _ in _molecule_groups(query))
    candidate_counts = Counter(label for label, _ in _molecule_groups(candidate))
    labels = sorted(set(query_counts) | set(candidate_counts))
    shared = [{"label": label, "count": min(query_counts[label], candidate_counts[label])}
              for label in labels if query_counts[label] and candidate_counts[label]]
    only_query = [{"label": label, "count": query_counts[label] - candidate_counts[label]}
                  for label in labels if query_counts[label] > candidate_counts[label]]
    only_candidate = [{"label": label, "count": candidate_counts[label] - query_counts[label]}
                      for label in labels if candidate_counts[label] > query_counts[label]]
    return {"shared": shared, "only_query": only_query, "only_candidate": only_candidate}


def _molecule_groups(mol):
    """Return named functional groups and rings as atom-index sets."""
    functional = []
    for priority, (label, smarts) in enumerate(_FUNCTIONAL_GROUPS):
        pattern = Chem.MolFromSmarts(smarts)
        if pattern is None:
            continue
        for match in mol.GetSubstructMatches(pattern, uniquify=True):
            if len(match) > 1 or label in ("Azole ring nitrogen", "Halogen (F/Cl/Br/I)"):
                functional.append((priority, label, frozenset(match)))
    groups = set()
    for priority, label, atoms in functional:
        shadowed = any(
            other_priority < priority
            and atoms.issubset(other_atoms)
            and (atoms != other_atoms or other_label != label)
            for other_priority, other_label, other_atoms in functional
        )
        if not shadowed:
            groups.add((label, atoms))
    for ring in mol.GetRingInfo().AtomRings():
        atoms = frozenset(ring)
        if atoms:
            ring_type = "Aromatic" if all(mol.GetAtomWithIdx(i).GetIsAromatic() for i in atoms) else "Non-aromatic"
            label = f"{ring_type} ring ({len(atoms)} atoms)"
            groups.add((label, atoms))
    return groups


def compare_molecular_features(query, candidate):
    """Summarize shared and side-specific functional groups and ring types."""
    query_counts = Counter(label for label, _ in _molecule_groups(query))
    candidate_counts = Counter(label for label, _ in _molecule_groups(candidate))
    labels = sorted(set(query_counts) | set(candidate_counts))
    shared = [{"label": label, "count": min(query_counts[label], candidate_counts[label])}
              for label in labels if query_counts[label] and candidate_counts[label]]
    only_query = [{"label": label, "count": query_counts[label] - candidate_counts[label]}
                  for label in labels if query_counts[label] > candidate_counts[label]]
    only_candidate = [{"label": label, "count": candidate_counts[label] - query_counts[label]}
                      for label in labels if candidate_counts[label] > query_counts[label]]
    return {"shared": shared, "only_query": only_query, "only_candidate": only_candidate}


def _molecule_groups(mol):
    """Return named functional groups and rings as atom-index sets."""
    functional = []
    for priority, (label, smarts) in enumerate(_FUNCTIONAL_GROUPS):
        pattern = Chem.MolFromSmarts(smarts)
        if pattern is None:
            continue
        for match in mol.GetSubstructMatches(pattern, uniquify=True):
            if len(match) > 1 or label in ("Azole ring nitrogen", "Halogen (F/Cl/Br/I)"):
                functional.append((priority, label, frozenset(match)))
    groups = set()
    for priority, label, atoms in functional:
        shadowed = any(
            other_priority < priority
            and atoms.issubset(other_atoms)
            and (atoms != other_atoms or other_label != label)
            for other_priority, other_label, other_atoms in functional
        )
        if not shadowed:
            groups.add((label, atoms))
    for ring in mol.GetRingInfo().AtomRings():
        atoms = frozenset(ring)
        if atoms:
            ring_type = "Aromatic" if all(mol.GetAtomWithIdx(i).GetIsAromatic() for i in atoms) else "Non-aromatic"
            label = f"{ring_type} ring ({len(atoms)} atoms)"
            groups.add((label, atoms))
    return groups


def compare_molecular_features(query, candidate):
    """Summarize shared and side-specific functional groups and ring types."""
    query_counts = Counter(label for label, _ in _molecule_groups(query))
    candidate_counts = Counter(label for label, _ in _molecule_groups(candidate))
    labels = sorted(set(query_counts) | set(candidate_counts))
    shared = [{"label": label, "count": min(query_counts[label], candidate_counts[label])}
              for label in labels if query_counts[label] and candidate_counts[label]]
    only_query = [{"label": label, "count": query_counts[label] - candidate_counts[label]}
                  for label in labels if query_counts[label] > candidate_counts[label]]
    only_candidate = [{"label": label, "count": candidate_counts[label] - query_counts[label]}
                      for label in labels if candidate_counts[label] > query_counts[label]]
    return {"shared": shared, "only_query": only_query, "only_candidate": only_candidate}


def pair_svgs(query, candidate):
    """Highlight a deterministic maximum common substructure with whole groups."""
    from rdkit.Chem import rdFMCS
    from rdkit.Chem.Draw import rdMolDraw2D
    common_query, common_candidate = set(), set()
    query_bonds, candidate_bonds = set(), set()
    qgroups, cgroups = _molecule_groups(query), _molecule_groups(candidate)
    qgroup_lookup, cgroup_lookup = set(qgroups), set(cgroups)
    try:
        try:
            mcs = rdFMCS.FindMCS([query, candidate], timeout=2,
                atomCompare=rdFMCS.AtomCompare.CompareElements,
                bondCompare=rdFMCS.BondCompare.CompareOrderExact,
                ringMatchesRingOnly=True, completeRingsOnly=True, matchValences=True)
        except (AttributeError, TypeError, ValueError):
            mcs = rdFMCS.FindMCS([query, candidate], timeout=2,
                atomCompare=rdFMCS.AtomCompare.CompareElements,
                bondCompare=rdFMCS.BondCompare.CompareOrder,
                ringMatchesRingOnly=True, completeRingsOnly=True, matchValences=True)
        pattern = Chem.MolFromSmarts(mcs.smartsString) if mcs.smartsString else None
        if pattern is not None:
            qmatches = query.GetSubstructMatches(pattern, uniquify=False, maxMatches=64)
            cmatches = candidate.GetSubstructMatches(pattern, uniquify=False, maxMatches=64)
            qrank = list(Chem.CanonicalRankAtoms(query, breakTies=True))
            crank = list(Chem.CanonicalRankAtoms(candidate, breakTies=True))
            qcanonical = Chem.MolToSmiles(query, canonical=True, isomericSmiles=True)
            ccanonical = Chem.MolToSmiles(candidate, canonical=True, isomericSmiles=True)
            best_key, best_pair, best_qgroups, best_cgroups = None, None, set(), set()
            for qmatch in qmatches:
                qset = set(qmatch)
                for cmatch in cmatches:
                    mapping = dict(zip(qmatch, cmatch))
                    shared_q, shared_c = set(), set()
                    for label, atoms in qgroup_lookup:
                        if atoms.issubset(qset):
                            mapped = frozenset(mapping[i] for i in atoms)
                            if (label, mapped) in cgroup_lookup:
                                shared_q.update(atoms)
                                shared_c.update(mapped)
                    if qcanonical <= ccanonical:
                        canonical_key = tuple(sorted((qrank[qi], crank[ci]) for qi, ci in zip(qmatch, cmatch)))
                    else:
                        canonical_key = tuple(sorted((crank[ci], qrank[qi]) for qi, ci in zip(qmatch, cmatch)))
                    key = (-len(shared_q), -len(shared_c), canonical_key)
                    if best_key is None or key < best_key:
                        best_key, best_pair = key, (qmatch, cmatch)
                        best_qgroups, best_cgroups = shared_q, shared_c
            if best_pair:
                qmatch, cmatch = best_pair
                q_group_atoms = set().union(*(set(atoms) for _, atoms in qgroups)) if qgroups else set()
                c_group_atoms = set().union(*(set(atoms) for _, atoms in cgroups)) if cgroups else set()
                common_query = (set(qmatch) - q_group_atoms) | best_qgroups
                common_candidate = (set(cmatch) - c_group_atoms) | best_cgroups
                for bond in pattern.GetBonds():
                    q1, q2 = qmatch[bond.GetBeginAtomIdx()], qmatch[bond.GetEndAtomIdx()]
                    c1, c2 = cmatch[bond.GetBeginAtomIdx()], cmatch[bond.GetEndAtomIdx()]
                    qbond = query.GetBondBetweenAtoms(q1, q2)
                    cbond = candidate.GetBondBetweenAtoms(c1, c2)
                    if qbond is not None and q1 in common_query and q2 in common_query:
                        query_bonds.add(qbond.GetIdx())
                    if cbond is not None and c1 in common_candidate and c2 in common_candidate:
                        candidate_bonds.add(cbond.GetIdx())
    except Exception:
        common_query.clear()
        common_candidate.clear()
        query_bonds.clear()
        candidate_bonds.clear()

    def draw(mol, shared_atoms, shared_bonds):
        atom_colors = {i: ((0.16, 0.75, 0.45) if i in shared_atoms else (0.93, 0.25, 0.35))
                       for i in range(mol.GetNumAtoms())}
        bond_colors = {b.GetIdx(): ((0.16, 0.75, 0.45) if b.GetIdx() in shared_bonds else (0.93, 0.25, 0.35))
                       for b in mol.GetBonds()}
        drawer = rdMolDraw2D.MolDraw2DSVG(360, 220)
        options = drawer.drawOptions()
        options.padding = 0.12
        options.clearBackground = False
        drawer.DrawMolecule(mol, highlightAtoms=list(atom_colors), highlightBonds=list(bond_colors),
            highlightAtomColors=atom_colors, highlightBondColors=bond_colors)
        drawer.FinishDrawing()
        svg = drawer.GetDrawingText()
        svg = re.sub(r"<rect\b(?=[^>]*(?:fill\s*:\s*#(?:fff|ffffff)\b|fill=['\"]#(?:fff|ffffff)['\"]))[^>]*(?:/>|>\s*</rect\s*>)", "", svg, flags=re.I)
        return svg

    return draw(query, common_query, query_bonds), draw(candidate, common_candidate, candidate_bonds)


def _compute_mcs_feature_match(query, candidate):
    """Return one deterministic MCS mapping and the exact groups shared by it."""
    from rdkit.Chem import rdFMCS

    qgroups, cgroups = _molecule_groups(query), _molecule_groups(candidate)
    qgroup_lookup, cgroup_lookup = set(qgroups), set(cgroups)
    try:
        try:
            mcs = rdFMCS.FindMCS([query, candidate], timeout=2,
                atomCompare=rdFMCS.AtomCompare.CompareElements,
                bondCompare=rdFMCS.BondCompare.CompareOrderExact,
                ringMatchesRingOnly=True, completeRingsOnly=True, matchValences=True)
        except (AttributeError, TypeError, ValueError):
            mcs = rdFMCS.FindMCS([query, candidate], timeout=2,
                atomCompare=rdFMCS.AtomCompare.CompareElements,
                bondCompare=rdFMCS.BondCompare.CompareOrder,
                ringMatchesRingOnly=True, completeRingsOnly=True, matchValences=True)
        pattern = Chem.MolFromSmarts(mcs.smartsString) if mcs.smartsString else None
        if pattern is None:
            return None
        qmatches = query.GetSubstructMatches(pattern, uniquify=False, maxMatches=64)
        cmatches = candidate.GetSubstructMatches(pattern, uniquify=False, maxMatches=64)
        qrank = list(Chem.CanonicalRankAtoms(query, breakTies=True))
        crank = list(Chem.CanonicalRankAtoms(candidate, breakTies=True))
        qcanonical = Chem.MolToSmiles(query, canonical=True, isomericSmiles=True)
        ccanonical = Chem.MolToSmiles(candidate, canonical=True, isomericSmiles=True)
        best_key, best = None, None
        for qmatch in qmatches:
            qset = set(qmatch)
            for cmatch in cmatches:
                mapping = dict(zip(qmatch, cmatch))
                matched_groups = []
                for label, atoms in qgroup_lookup:
                    if atoms.issubset(qset):
                        mapped = frozenset(mapping[index] for index in atoms)
                        if (label, mapped) in cgroup_lookup:
                            matched_groups.append((label, atoms, mapped))
                shared_q = set().union(*(set(atoms) for _, atoms, _ in matched_groups)) if matched_groups else set()
                shared_c = set().union(*(set(atoms) for _, _, atoms in matched_groups)) if matched_groups else set()
                if qcanonical <= ccanonical:
                    canonical_key = tuple(sorted((qrank[qi], crank[ci]) for qi, ci in zip(qmatch, cmatch)))
                else:
                    canonical_key = tuple(sorted((crank[ci], qrank[qi]) for qi, ci in zip(qmatch, cmatch)))
                key = (-len(shared_q), -len(shared_c), canonical_key)
                if best_key is None or key < best_key:
                    best_key = key
                    best = (pattern, qmatch, cmatch, matched_groups)
        if best is None:
            return None
        pattern, qmatch, cmatch, matched_groups = best
        q_group_atoms = set().union(*(set(atoms) for _, atoms in qgroups)) if qgroups else set()
        c_group_atoms = set().union(*(set(atoms) for _, atoms in cgroups)) if cgroups else set()
        matched_q = set().union(*(set(atoms) for _, atoms, _ in matched_groups)) if matched_groups else set()
        matched_c = set().union(*(set(atoms) for _, _, atoms in matched_groups)) if matched_groups else set()
        shared_query = (set(qmatch) - q_group_atoms) | matched_q
        shared_candidate = (set(cmatch) - c_group_atoms) | matched_c
        query_bonds, candidate_bonds = set(), set()
        for bond in pattern.GetBonds():
            q1, q2 = qmatch[bond.GetBeginAtomIdx()], qmatch[bond.GetEndAtomIdx()]
            c1, c2 = cmatch[bond.GetBeginAtomIdx()], cmatch[bond.GetEndAtomIdx()]
            qbond = query.GetBondBetweenAtoms(q1, q2)
            cbond = candidate.GetBondBetweenAtoms(c1, c2)
            if qbond is not None and q1 in shared_query and q2 in shared_query:
                query_bonds.add(qbond.GetIdx())
            if cbond is not None and c1 in shared_candidate and c2 in shared_candidate:
                candidate_bonds.add(cbond.GetIdx())
        return {
            "shared_query": shared_query, "shared_candidate": shared_candidate,
            "query_bonds": query_bonds, "candidate_bonds": candidate_bonds,
            "query_groups": qgroups, "candidate_groups": cgroups,
            "matched_groups": matched_groups,
        }
    except Exception:
        return None


_MCS_FEATURE_CACHE = {}


def _mcs_feature_match(query, candidate):
    """Reuse a pair's mapping for both the SVG and its side-specific feature list."""
    try:
        key = (
            Chem.MolToSmiles(query, canonical=True, isomericSmiles=True),
            Chem.MolToSmiles(candidate, canonical=True, isomericSmiles=True),
        )
    except Exception:
        return _compute_mcs_feature_match(query, candidate)
    if key in _MCS_FEATURE_CACHE:
        return _MCS_FEATURE_CACHE[key]
    match = _compute_mcs_feature_match(query, candidate)
    if len(_MCS_FEATURE_CACHE) >= 256:
        _MCS_FEATURE_CACHE.pop(next(iter(_MCS_FEATURE_CACHE)))
    _MCS_FEATURE_CACHE[key] = match
    return match


def _unmatched_fragments(mol, shared_atoms, covered_atoms):
    """Describe red fragments that are not already represented by a named group."""
    remaining = set(range(mol.GetNumAtoms())) - set(shared_atoms) - set(covered_atoms)
    fragments = []
    while remaining:
        stack = [min(remaining)]
        component = set()
        while stack:
            index = stack.pop()
            if index in component or index not in remaining:
                continue
            component.add(index)
            stack.extend(neighbor.GetIdx() for neighbor in mol.GetAtomWithIdx(index).GetNeighbors()
                         if neighbor.GetIdx() in remaining)
        remaining.difference_update(component)
        smiles = Chem.MolFragmentToSmiles(mol, atomsToUse=sorted(component), canonical=True, isomericSmiles=True)
        attachments = sorted({
            neighbor.GetSymbol()
            for index in component
            for neighbor in mol.GetAtomWithIdx(index).GetNeighbors()
            if neighbor.GetIdx() in shared_atoms
        })
        attached = f"; attached to {'/'.join(attachments)}" if attachments else ""
        if smiles:
            fragments.append({"label": f"Unmatched fragment: {smiles} ({len(component)} atoms{attached})", "count": 1})
    return fragments


def _feature_rows(groups):
    counts = Counter(label for label, _ in groups)
    return [{"label": label, "count": count} for label, count in sorted(counts.items())]


def compare_molecular_features(query, candidate):
    """List groups only when the MCS mapping confirms they are actually shared."""
    match = _mcs_feature_match(query, candidate)
    qgroups, cgroups = _molecule_groups(query), _molecule_groups(candidate)
    if match is None:
        shared_groups, only_qgroups, only_cgroups = [], [], []
        shared_query, shared_candidate = set(), set()
    else:
        matched = match["matched_groups"]
        shared_groups = [(label, qatoms) for label, qatoms, _ in matched]
        matched_q = {(label, atoms) for label, atoms, _ in matched}
        matched_c = {(label, atoms) for label, _, atoms in matched}
        only_qgroups = [group for group in qgroups if group not in matched_q]
        only_cgroups = [group for group in cgroups if group not in matched_c]
        shared_query, shared_candidate = match["shared_query"], match["shared_candidate"]
    only_q_atoms = set().union(*(set(atoms) for _, atoms in only_qgroups)) if only_qgroups else set()
    only_c_atoms = set().union(*(set(atoms) for _, atoms in only_cgroups)) if only_cgroups else set()
    return {
        "shared": _feature_rows(shared_groups),
        "only_query": _feature_rows(only_qgroups) + _unmatched_fragments(query, shared_query, only_q_atoms),
        "only_candidate": _feature_rows(only_cgroups) + _unmatched_fragments(candidate, shared_candidate, only_c_atoms),
    }


def pair_svgs(query, candidate):
    """Highlight the same shared groups and unmatched atoms used by the text summary."""
    from rdkit.Chem.Draw import rdMolDraw2D

    match = _mcs_feature_match(query, candidate) or {}

    def draw(mol, shared_atoms, shared_bonds):
        atom_colors = {i: ((0.16, 0.75, 0.45) if i in shared_atoms else (0.93, 0.25, 0.35))
                       for i in range(mol.GetNumAtoms())}
        bond_colors = {b.GetIdx(): ((0.16, 0.75, 0.45) if b.GetIdx() in shared_bonds else (0.93, 0.25, 0.35))
                       for b in mol.GetBonds()}
        drawer = rdMolDraw2D.MolDraw2DSVG(360, 220)
        options = drawer.drawOptions()
        options.padding = 0.12
        options.clearBackground = False
        drawer.DrawMolecule(mol, highlightAtoms=list(atom_colors), highlightBonds=list(bond_colors),
            highlightAtomColors=atom_colors, highlightBondColors=bond_colors)
        drawer.FinishDrawing()
        svg = drawer.GetDrawingText()
        return re.sub(r"<rect\b(?=[^>]*(?:fill\s*:\s*#(?:fff|ffffff)\b|fill=['\"]#(?:fff|ffffff)['\"]))[^>]*(?:/>|>\s*</rect\s*>)", "", svg, flags=re.I)

    return (
        draw(query, match.get("shared_query", set()), match.get("query_bonds", set())),
        draw(candidate, match.get("shared_candidate", set()), match.get("candidate_bonds", set())),
    )
