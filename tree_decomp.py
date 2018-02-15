import itertools
from rdkit import Chem
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
  return all_clusters


def create_substructure(m, atom_ids):
  emol = EditableMol(m)
  to_remove_ids = [x.GetIdx() for x in m.GetAtoms()]
  for atom_id in atom_ids:
    to_remove_ids.remove(atom_id)
  to_remove_ids.reverse()
  for atom_id in to_remove_ids:
    emol.RemoveAtom(atom_id)
  return emol.GetMol()


def main():
  lines = open('data/250k_rndm_zinc_drugs_clean.smi').readlines()
  mols = [Chem.MolFromSmiles(x) for x in lines]
  substructure_smiles = set()
  for mol in mols:
    clusters = get_cluster_atoms(mol)
    for cluster in clusters:
      fragment = create_substructure(mol, cluster)
      substructure_smiles.add(Chem.MolToSmiles(fragment))
  print(len(substructure_smiles))
  with open('data/fragments.txt', 'w') as fout:
    fout.write(json.dumps(list(substructure_smiles)))


if __name__ == "__main__":
  main()
