"""SA fingerprint attribution and atom-aligned illustrative 3D conformer."""
import math
from functools import lru_cache
from rdkit import Chem
from rdkit.Chem import AllChem, rdMolDescriptors
from .synth import _load_sascorer

VERSION = 'sa3d-1'

@lru_cache(maxsize=32)
def explain(smiles):
    mol = Chem.MolFromSmiles(smiles.strip())
    if mol is None or not mol.GetNumAtoms():
        raise ValueError('Enter a valid, non-empty SMILES.')
    mol = Chem.RemoveHs(mol)
    if mol.GetNumAtoms() > 300:
        raise ValueError('The interactive explanation supports up to 300 heavy atoms.')
    sas = _load_sascorer()
    if sas is None:
        raise ValueError('The RDKit SA scorer is unavailable in this runtime.')
    actual = float(sas.calculateScore(mol))
    # Legacy sparse Morgan interface supplies the same unfurled IDs and bitInfo
    # across older and newer RDKit releases. Verify against installed scorer below.
    info = {}
    fp = AllChem.GetMorganFingerprint(mol, 2, bitInfo=info)
    counts = fp.GetNonzeroElements()
    if hasattr(sas, 'mfpgen'):
        expected = sas.mfpgen.GetSparseCountFingerprint(mol).GetNonzeroElements()
        if dict(counts) != dict(expected):
            raise ValueError('Installed scorer uses different Morgan settings; explanation stopped to avoid incorrect atom highlights.')
    total = sum(counts.values())
    rows = []
    for bit,count in counts.items():
        weight = float(sas._fscores.get(bit, -4.0))
        environments = []
        for centre,radius in info[bit]:
            bonds = list(Chem.FindAtomEnvironmentOfRadiusN(mol, int(radius), int(centre))) if radius else []
            atoms = {int(centre)}
            for bi in bonds:
                bond = mol.GetBondWithIdx(bi)
                atoms.update([bond.GetBeginAtomIdx(),bond.GetEndAtomIdx()])
            centre_atom = mol.GetAtomWithIdx(int(centre))
            if centre_atom.GetIsAromatic():
                centre_label = f"aromatic {centre_atom.GetSymbol()}"
            elif centre_atom.GetAtomicNum() == 6 and any(
                    b.GetBondType() == Chem.BondType.DOUBLE and b.GetOtherAtom(centre_atom).GetAtomicNum() == 8
                    for b in centre_atom.GetBonds()):
                centre_label = "carbonyl carbon"
            elif centre_atom.GetAtomicNum() == 8 and any(
                    b.GetBondType() == Chem.BondType.DOUBLE and b.GetOtherAtom(centre_atom).GetAtomicNum() == 6
                    for b in centre_atom.GetBonds()):
                centre_label = "carbonyl oxygen"
            else:
                centre_label = centre_atom.GetSymbol()
            charge_label = f"charge {centre_atom.GetFormalCharge():+d}" if centre_atom.GetFormalCharge() else "neutral"
            centre_label = (f"{centre_label} (degree {centre_atom.GetDegree()}, valence "
                            f"{centre_atom.GetTotalValence()}, H {centre_atom.GetTotalNumHs()}, {charge_label})")
            fragment_label = f"{centre_label} environment · radius {int(radius)}"
            fragment_smiles = Chem.MolFragmentToSmiles(
                mol, atomsToUse=sorted(atoms), bondsToUse=bonds,
                canonical=True, isomericSmiles=True)
            environments.append({'centre':int(centre),'radius':int(radius),
                                 'atoms':sorted(atoms),'bonds':bonds,
                                 'fragment_smiles':fragment_smiles,
                                 'fragment_label':fragment_label})
        rows.append({'id':str(bit),'count':int(count),'matched':bit in sas._fscores,
                     'weight':weight,'contribution':count*weight/total,'environments':environments})
    rows.sort(key=lambda x:(-abs(x['contribution']),x['id']))
    n = mol.GetNumAtoms()
    penalties = {'Size':n**1.005-n,
        'Stereocentres':math.log10(len(Chem.FindMolChiralCenters(mol,includeUnassigned=True))+1),
        'Spiro atoms':math.log10(rdMolDescriptors.CalcNumSpiroAtoms(mol)+1),
        'Bridgehead atoms':math.log10(rdMolDescriptors.CalcNumBridgeheadAtoms(mol)+1),
        'Macrocycles':math.log10(2) if any(len(r)>8 for r in mol.GetRingInfo().AtomRings()) else 0.0}
    fragment = sum(x['contribution'] for x in rows)
    penalty = sum(penalties.values())
    correction = math.log(n/len(rows))*0.5 if n>len(rows) else 0.0
    raw = fragment-penalty+correction
    mapped = 11-(raw+5)*9/6.5
    reconstructed = max(1,min(10,8+math.log(mapped-8) if mapped>8 else mapped))
    if abs(reconstructed-actual)>1e-7:
        raise ValueError('Explanation does not reproduce this installed scorer version. SA calculation remains unchanged.')
    # Add/remove Hs without canonicalising/reordering the heavy atoms.
    display_mol = Chem.AddHs(Chem.Mol(mol))
    params = AllChem.ETKDGv3() if hasattr(AllChem,'ETKDGv3') else AllChem.ETKDGv2()
    params.randomSeed = 42
    params.maxIterations = 200
    if hasattr(params,'timeout'):params.timeout = 15
    conformer_ok = AllChem.EmbedMolecule(display_mol,params)==0
    molblock = None
    if conformer_ok:
        display_mol = Chem.RemoveHs(display_mol)
        # Atom metadata is compared before serializing; SDF preserves this order.
        if [a.GetAtomicNum() for a in mol.GetAtoms()] != [a.GetAtomicNum() for a in display_mol.GetAtoms()]:
            raise ValueError('3D atom mapping validation failed.')
        molblock = Chem.MolToMolBlock(display_mol)
    sa = round(actual,2)
    return {'version':VERSION,'smiles':Chem.MolToSmiles(mol),'sa':sa,
            'score':int(round(round((10-sa)/9,4)*100)),
            'rows':rows,'unique':len(rows),'occurrences':total,
            'matched_unique':sum(x['matched'] for x in rows),
            'matched_occurrences':sum(x['count'] for x in rows if x['matched']),
            'fragment_mean':fragment,'penalties':penalties,'penalty_total':penalty,
            'density_correction':correction,'raw':raw,'mapped':mapped,'verified':True,
            'molblock':molblock,'atom_count':n,'elements':[a.GetSymbol() for a in mol.GetAtoms()],
            'viewer_note':'Generated conformer for illustration, not a docked pose. SA uses 2D structure.' if conformer_ok else '3D embedding failed. The verified fragment table and SA calculation remain available.'}
