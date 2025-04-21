from ring_search import ring_counting
from graph_creation import max_ring_adaptation, create_graph_from_mol

import numpy as np
import matplotlib.pyplot as plt
import argparse

from ase.io import read

parser = argparse.ArgumentParser()
parser.add_argument("--xyz")
parser.add_argument("--pbc", default=True)
parser.add_argument("--ring", default=20)
args = parser.parse_args()

def main():
    system = read(args.xyz)

    graph_system = max_ring_adaptation(system, periodic=args.pbc, max_ring=int(args.ring))
    molecular_graph = create_graph_from_mol(molecule=graph_system)

    atoms = [atom.index for atom in system if atom.symbol == "Si"]
    rings, ring_amount = ring_counting(T_atoms=atoms, silica_graph=molecular_graph, max_ring=int(args.ring))
    print(np.array(list(ring_amount.values()))/len(atoms))
    # plt.bar(ring_amount.keys(), np.array(list(ring_amount.values()))/len(atoms))
    # plt.show()

if __name__ == "__main__":
    main()
