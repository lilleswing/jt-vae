import itertools
from rdkit import Chem
import copy
from rdkit.Chem import EditableMol
import networkx as nx
import json

flatten = lambda l: [item for sublist in l for item in sublist]


def get_cluster_atoms(m):
  """
  NOTE(LESWING) v0 calculation is needed for maximum spanning tree calculation
  :param m:
  :return:
  """
  r_info = m.GetRingInfo()
  bond_rings = set(flatten(r_info.BondRings()))
  all_bonds = [x.GetIdx() for x in m.GetBonds()]
  non_ring_bonds = set(all_bonds) - bond_rings
  v1 = set([tuple(sorted(
    [m.GetBonds()[x].GetBeginAtom().GetIdx(), m.GetBonds()[x].GetEndAtom().GetIdx()]))
    for x in non_ring_bonds])

  v2 = set()
  v2.update([tuple(sorted((x))) for x in r_info.AtomRings()])

  # Merge Rings
  to_merge = set()
  for r1, r2 in itertools.product(v2, repeat=2):
    if r1 >= r2:
      continue
    intersection = set(r1).intersection(set(r2))
    if len(intersection) >= 3:
      to_merge.add((r1, r2))

  g = nx.Graph()
  for f, t in to_merge:
    g.add_edge(f, t)
  graphs = list(nx.connected_component_subgraphs(g))
  to_merge = [list(x.nodes()) for x in graphs]

  for merge_keys in to_merge:
    s1 = set()
    for merge_key in merge_keys:
      v2.remove(merge_key)
      s1.update(merge_key)
    s1 = tuple(sorted(list(s1)))
    v2.add(s1)

  all_clusters = set()
  all_clusters.update(v1)
  all_clusters.update(v2)
  v0 = set()
  for atom in m.GetAtoms():
    atom_id = atom.GetIdx()
    count = sum([
      atom_id in x for x in all_clusters
    ])
    if count >= 3:
      v0.add(atom_id)
  v1 = [('non_ring', x) for x in v1]
  v2 = [('ring', x) for x in v2]
  all_clusters = set()
  all_clusters.update(v1)
  all_clusters.update(v2)
  return all_clusters


def bonds_connected_to_atom(mol, atom_idx):
  atoms_used = set([])
  queue = [atom_idx]
  bonds_used = set()
  while len(queue) > 0:
    node, queue = queue[0], queue[1:]
    if node in atoms_used:
      continue
    atoms_used.add(node)
    atom = mol.GetAtomWithIdx(node)
    bonds = atom.GetBonds()
    for bond in bonds:
      begin, end = bond.GetBeginAtom().GetIdx(), bond.GetEndAtom().GetIdx()
      bonds_used.add(bond.GetIdx())
      if begin not in atoms_used:
        queue.append(begin)
      if end not in atoms_used:
        queue.append(end)
  return list(bonds_used)


def frag_on_bonds(m, atom_ids):
  breaking_bonds = list()
  for bond in m.GetBonds():
    start, end = bond.GetBeginAtom().GetIdx(), bond.GetEndAtom().GetIdx()
    if start in atom_ids and end not in atom_ids:
      breaking_bonds.append(bond.GetIdx())
    if end in atom_ids and start not in atom_ids:
      breaking_bonds.append(bond.GetIdx())
  my_copy = copy.deepcopy(m)
  return Chem.FragmentOnBonds(my_copy, breaking_bonds)


def replace_dummy_with_h(m):
  emol = EditableMol(m)
  for atom in m.GetAtoms():
    if atom.GetAtomicNum() == 1:
      emol.ReplaceAtom(atom.GetIdx(), Chem.Atom(0))
  m = emol.GetMol()
  # for atom in m.GetAtoms():
  #   if atom.GetAtomicNum() == 0:
  #     for bond in atom.GetBonds():
  #       bond.SetBondType(Chem.BondType.UNSPECIFIED)
    # atom.SetIsAromatic(False)
  return m


def update_aromatic(m, cluster):
  if cluster[0] == 'non_ring':
    for atom in m.GetAtoms():
      atom.SetIsAromatic(False)
  print(Chem.MolToSmiles(m), cluster)
  Chem.Kekulize(m, clearAromaticFlags=True)
  return m


def fast_break(mol, atom_ids):
  atom_ids = atom_ids[1]
  bonds_to_keep = []
  for bond in mol.GetBonds():
    if bond.GetBeginAtom().GetIdx() in atom_ids and bond.GetEndAtom().GetIdx() in atom_ids:
      bonds_to_keep.append(bond.GetIdx())
  return Chem.PathToSubmol(mol, bonds_to_keep)


def create_substructure(mol, cluster):
  try:
    split_mol = frag_on_bonds(mol, cluster[1])
    bonds_to_keep = bonds_connected_to_atom(split_mol, list(cluster)[1][0])
    frag = Chem.PathToSubmol(split_mol, bonds_to_keep)
    frag.UpdatePropertyCache()
    return replace_dummy_with_h(frag)
  except Exception as e:
    print(Chem.MolToSmiles(mol))
    # Giant single "cluster" things
    # C1CSCCSCCS1
    # C1CSCCSCCCSCCSC1
    return None


def main():
  lines = open('data/250k_rndm_zinc_drugs_clean.smi').readlines()
  mols = [x for x in lines]
  substructure_smiles = set()
  for mol in mols:
    mol = Chem.MolFromSmiles(mol)
    Chem.Kekulize(mol)
    Chem.AddHs(mol)
    clusters = get_cluster_atoms(mol)
    for cluster in clusters:
      fragment = fast_break(mol, cluster)
      if fragment is None:
        continue
      substructure_smiles.add(Chem.MolToSmiles(fragment, kekuleSmiles=True))
  print(len(substructure_smiles))
  with open('data/fragments.txt', 'w') as fout:
    fout.write(json.dumps(list(substructure_smiles)))


if __name__ == "__main__":
  main()
