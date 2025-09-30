from rdkit import Chem
from rdkit.Chem import rdchem
import hashlib

try:
    from rdkit.Chem import rdChemReactions as _rxnmod
except Exception:
    from rdkit.Chem import AllChem as _rxnmod

def _rxn_from_smarts(s):
    try:
        return _rxnmod.ReactionFromSmarts(str(s), useSmiles=False)
    except Exception:
        return None

def _bond_set(mols):
    """Return set of (mapA,mapB,bondType,isAromatic) for all bonds in a list of mols.
       Only includes bonds whose atoms both have map numbers."""
    out = set()
    for m in mols:
        for b in m.GetBonds():
            a = b.GetBeginAtom(); z = b.GetEndAtom()
            ma = a.GetAtomMapNum(); mz = z.GetAtomMapNum()
            if ma and mz:
                a1, a2 = (ma, mz) if ma < mz else (mz, ma)
                out.add((a1, a2, int(b.GetBondType()), b.GetIsAromatic()))
    return out

def _reacting_maps(rxn):
    """Map numbers of atoms participating in any bond change."""
    r_bonds = _bond_set(list(rxn.GetReactants()))
    p_bonds = _bond_set(list(rxn.GetProducts()))
    changed = (r_bonds ^ p_bonds)
    maps = set()
    for a1, a2, *_ in changed:
        maps.add(a1); maps.add(a2)
    return maps, r_bonds, p_bonds

def _first_shell_signature(prod_mols, center_maps):
    """Signature = tuple of atoms (with features) and bonds in subgraph:
       center atoms + neighbors (radius 1) on product side."""
    idx = {}
    for m in prod_mols:
        for a in m.GetAtoms():
            amap = a.GetAtomMapNum()
            if amap:
                idx[amap] = (m, a.GetIdx())

    atom_keys = []
    bond_keys = set()

    def atom_feat(a: rdchem.Atom):
        try:
            hyb = int(a.GetHybridization())
        except Exception:
            hyb = -1
        return (
            a.GetSymbol(),
            a.GetIsAromatic(),
            a.IsInRing(),
            a.GetFormalCharge(),
            hyb,
        )

    queue = set(center_maps)
    for cmap in list(center_maps):
        m, ai = idx.get(cmap, (None, None))
        if m is None:
            continue
        a = m.GetAtomWithIdx(ai)
        for nb in a.GetNeighbors():
            if nb.GetAtomMapNum():
                queue.add(nb.GetAtomMapNum())

    for amap in queue:
        if amap in idx:
            m, ai = idx[amap]
            a = m.GetAtomWithIdx(ai)
            atom_keys.append((amap, atom_feat(a)))

    amap_set = set(queue)
    for m in prod_mols:
        for b in m.GetBonds():
            a = b.GetBeginAtom(); z = b.GetEndAtom()
            ma = a.GetAtomMapNum(); mz = z.GetAtomMapNum()
            if ma and mz and (ma in amap_set) and (mz in amap_set):
                a1, a2 = (ma, mz) if ma < mz else (mz, ma)
                bond_keys.add((a1, a2, int(b.GetBondType()), b.GetIsAromatic()))

    atoms_canon = tuple((feat) for _, feat in sorted(atom_keys, key=lambda x: x[0]))
    bonds_canon = tuple(sorted([(t, aro) for _, _, t, aro in bond_keys]))

    return atoms_canon, bonds_canon

def compute_center_group_key(template_smarts: str) -> str:
    """
    group key: exact reacting center + first shell on product side AND bond-change multiset.
    """
    rxn = _rxn_from_smarts(template_smarts)
    if rxn is None:
        return hashlib.sha256(("RAW|" + template_smarts).encode()).hexdigest()

    center_maps, r_bonds, p_bonds = _reacting_maps(rxn)

    formed = tuple(sorted(p_bonds - r_bonds))
    broken = tuple(sorted(r_bonds - p_bonds))
    changed = tuple(sorted(r_bonds & p_bonds))

    atoms_sig, bonds_sig = _first_shell_signature(list(rxn.GetProducts()), center_maps)

    payload = ("L1|", atoms_sig, bonds_sig, formed, broken, changed)
    raw = repr(payload).encode()
    return hashlib.sha256(raw).hexdigest()
