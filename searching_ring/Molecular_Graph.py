from ase import Atoms
from ase.neighborlist import NeighborList, natural_cutoffs
from ase.build.supercells import make_supercell

import networkx as nx
import numpy as np

from itertools import combinations

class Molecule_To_Graph():
    '''
    Class to create a networkX Graph object from an ASE atoms object.
    It is useful to apply graph theory to molecules.
    ----------------------------------------------------------------------------------
    Attributes:
            structure: ASE atoms object of the molecule
            periodic: 3x1 np.array describing the structure periodicity [True, True, True] for a 3D periodic
            box_size: 6x1 np.array of the molecule unit cell [a b c $\alpha$ $\beta$ $\gamma$]
            silica: bool specify if it is a silica gramework (allow to reduce the graph by half)
            max_ring: largest ring size to look for
            graph: networkX graph object of the molecule
            nodes: list[int] selected nodes from the graph (e.g. list of nodes to perform primitive ring search)
    '''
    def __init__(self, structure:Atoms, periodic:np.array=np.full(3, False), box_size:np.array=None, silica:bool=False, max_ring:int=20, nodes:list=None):
        self.structure = structure
        self.periodic = periodic
        self.box_size = box_size
        self.silica = silica
        self.max_ring = max_ring
        self.graph = nx.Graph()
        self.nodes = nodes if nodes is not None else []

    def max_ring_adaptation(self):
        '''
        Create a supercell from a periodic silica structure if one of the edges is smaller than the maximum ring size.
        The maximum ring size is computed as (max_ring_size/2)/3.4
        '''
        min_length = self.max_ring/2*1.6
        cell_extension = min_length//self.structure.cell.lengths()
        if  np.allclose(cell_extension, np.zeros(3)) is False:
            size_extension = np.diag(cell_extension + np.ones(3))
            self.structure = make_supercell(self.structure, P=size_extension, wrap=False)
            print(f"Structure too small to compute every {self.max_ring}-rings, must create a {np.diag(size_extension)} supercell.")

    def reduced_graph(self):
        '''
        Create a smaller graph from a graph describing a silica structure.
        Remove the oxygen atoms and use them to describing edges in the graph.
        '''
        si_index = [atom.index for atom in self.structure if atom.symbol == "Si"]
        edges = []
        silicon_graph = nx.Graph()
        silicon_graph.add_nodes_from(si_index)
        for oxygen in self.structure:
            if oxygen.symbol in ["O"]:
                neighbor_list = list(self.graph.neighbors(oxygen.index))
                if len(neighbor_list) > 1:
                    for pairs in combinations(neighbor_list, 2):
                        edges.append(pairs) 
        silicon_graph.add_edges_from(edges)
        self.graph = silicon_graph
        new_nodes = [index for index in self.nodes if self.structure[index].symbol == "Si"]
        self.nodes = new_nodes
    
    def create_graph_from_mol(self):
        '''
        Create a graph from an ASE atoms object
        '''
        self.structure.pbc[False, False, False]
        cutoffs = natural_cutoffs(self.structure)
        neighbor_list = NeighborList(cutoffs, self_interaction=False, bothways=True)
        neighbor_list.update(self.structure)
        matrix = neighbor_list.get_connectivity_matrix()
        self.graph = nx.from_scipy_sparse_array(matrix)
        if len(self.nodes) == 0:
            self.nodes = range(len(self.graph))

    def super_cell_periodicity(self):
        '''
        To perform primitive ring search on periodic system, we need to do it on a supercell to find rings at the frontiers
        '''
        self.structure.pbc = self.periodic
        if self.box_size is not None:
            self.structure.cell = np.array(self.box_size)
        self.max_ring_adaptation()
        matrix = np.diag(self.periodic)*np.full(3, 3)
        dimensionality = np.sum(np.array(self.periodic))
        replica_to_chose = 3**dimensionality//2
        self.nodes = [x + len(self.structure.get_atomic_numbers())*replica_to_chose for x in range(len(self.structure.get_atomic_numbers()))]
        self.structure = make_supercell(prim=self.structure, P=matrix)

    def full_routine(self):
        '''
        Main routine to create the graph form the ASE atoms.
        Knowing if the system is periodic or a silica structure.
        '''
        if np.any(self.periodic):
            self.super_cell_periodicity()
        self.create_graph_from_mol()
        if self.silica == True:
            self.reduced_graph() 
        return self.graph, self.nodes

