from ase import Atoms
from ase.neighborlist import NeighborList, natural_cutoffs
from ase.build.supercells import make_supercell

import networkx as nx
import numpy as np
import itertools

def create_graph_from_mol(molecule:Atoms):
    cutoffs = natural_cutoffs(molecule)
    neighbor_list = NeighborList(cutoffs, self_interaction=False, bothways=True)
    neighbor_list.update(molecule)
    matrix = neighbor_list.get_connectivity_matrix()
    return nx.from_scipy_sparse_array(matrix)

def max_ring_adaptation(molecule:Atoms, periodic:bool=True, max_ring:int=20):
    extended_system = molecule
    if periodic is False:
        molecule.set_pbc([False, False, False])
    elif periodic is True:
        min_length = max_ring/2*1.6
        cell_extension = min_length//molecule.cell.lengths()
        if  np.allclose(cell_extension, np.zeros(3)) is False:
            size_extension = np.diag(cell_extension + np.ones(3))
            extended_system = make_supercell(molecule, P=size_extension, wrap=False)
            print(f"Structure too small to compute every {max_ring}-rings, must create a {np.diag(size_extension)} supercell.")
    return extended_system

def reduced_graph(silica_graph:nx.Graph, molecule:Atoms):
    si_index = [atom.index for atom in molecule if atom.symbol == "Si"]
    edges = []
    silicon_graph = nx.Graph()
    silicon_graph.add_nodes_from(si_index)
    for oxygen in molecule:
        if oxygen.symbol == "O":
            neighbor_list = list(silica_graph.neighbors(oxygen.index))
            if len(neighbor_list) > 1:
                for pairs in itertools.combinations(neighbor_list, 2):
                    edges.append(pairs) 
    silicon_graph.add_edges_from(edges)
    return silicon_graph