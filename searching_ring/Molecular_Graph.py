from ase import Atoms
from ase.neighborlist import NeighborList, natural_cutoffs
from ase.build.supercells import make_supercell

import networkx as nx
import numpy as np

from itertools import combinations

class Molecule_To_Graph():
    def __init__(self, structure:Atoms, periodic:bool=False, box_size:np.array=None, silica:bool=False, max_ring:int=20, nodes:list=None):
        self.structure = structure
        self.periodic = periodic
        self.box_size = box_size
        self.silica = silica
        self.max_ring = max_ring
        self.graph = nx.Graph()
        self.nodes = nodes if nodes is not None else []

    def max_ring_adaptation(self):
        min_length = self.max_ring/2*1.6
        cell_extension = min_length//self.structure.cell.lengths()
        if  np.allclose(cell_extension, np.zeros(3)) is False:
            size_extension = np.diag(cell_extension + np.ones(3))
            self.structure = make_supercell(self.structure, P=size_extension, wrap=False)
            print(f"Structure too small to compute every {self.max_ring}-rings, must create a {np.diag(size_extension)} supercell.")

    def reduced_graph(self):
        si_index = [atom.index for atom in self.structure if atom.symbol == "Si"]
        edges = []
        silicon_graph = nx.Graph()
        silicon_graph.add_nodes_from(si_index)
        for oxygen in self.structure:
            if oxygen.symbol == "O":
                neighbor_list = list(self.graph.neighbors(oxygen.index))
                if len(neighbor_list) > 1:
                    for pairs in combinations(neighbor_list, 2):
                        edges.append(pairs) 
        silicon_graph.add_edges_from(edges)
        self.graph = silicon_graph
        self.nodes = self.graph.nodes
    
    def create_graph_from_mol(self):
        self.structure.pbc[False, False, False]
        cutoffs = natural_cutoffs(self.structure)
        neighbor_list = NeighborList(cutoffs, self_interaction=False, bothways=True)
        neighbor_list.update(self.structure)
        matrix = neighbor_list.get_connectivity_matrix()
        self.graph = nx.from_scipy_sparse_array(matrix)
        self.nodes = range(len(self.graph))

    def super_cell_periodicity(self):
        if self.box_size is not None:
            self.structure.cell = np.array(self.box_size)
            self.structure.pbc = [True, True, True]
        print(self.structure.cell, self.structure.pbc)
        self.max_ring_adaptation()
        matrix = np.array([[3, 0, 0], [0, 3, 0], [0, 0, 3]])
        self.structure = make_supercell(prim=self.structure, P=matrix)
        self.nodes = [x + len(self.structure.get_atomic_numbers())*13 for x in range(len(self.structure.get_atomic_numbers()))]

    def full_routine(self):
        if self.periodic == True:
            self.super_cell_periodicity()
        self.create_graph_from_mol()
        if self.silica == True:
            self.reduced_graph() 
        return self.graph, self.nodes

