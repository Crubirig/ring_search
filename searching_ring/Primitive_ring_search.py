import networkx as nx
import numpy as np
import matplotlib.pyplot as plt
import time

from itertools import combinations
from collections import deque

class Primitive_ring_search():
    '''
    Class to search primitive rings in a graph.
    -----------------------------------------------------------------------------------------------
    Attributes:
            graph: networkx Graph object
            max_size: largest ring size to look for
            src_nodes: list[int] of nodes for which to search for primitive rings
            distances: dict[dict[int]] containing the shortest path length between each pair of nodes
            prim_rings: dict[list] containing every primitive rings computed
    '''
    def __init__(self, graph:nx.Graph, max_size:int, distances:dict[dict]=None, src_nodes:list=[0], prim_rings:dict=None):
        self.graph = graph
        self.max_size = max_size
        self.src_nodes = src_nodes
        self.distances = distances if distances is not None else dict(
            nx.all_pairs_shortest_path_length(graph, cutoff=max_size//2 + 1))
        self.prim_rings = prim_rings if prim_rings is not None else {x:[] for x in range(2, self.get_max_size() + 2)}

    def get_graph(self):
        return self.graph

    def get_max_size(self):
        return self.max_size
    
    def get_shortest_path_lengths(self):
        return self.distances

    def find_all_shortest_paths(self, starting_node: int, final_node: int) -> list[list]:
        '''
        Find all shortest paths between two nodes in a graph, knowing the length of those paths. 
        This is much faster than when the length is not known.
        -------------------------------------------------------------------------------------------
        Inputs:
            graph: networkx Graph object
            starting_node: starting node index
            final_node: target node index
            shortest_path_length: dict[dict[int]] containing shortest paths length for each node pair

        Outputs:
            A list of every shortest path between the starting and target node
        '''

        queue = deque([[starting_node]])  # Using deque for better performance on queue operations
        paths = []
        if final_node not in self.get_shortest_path_lengths()[starting_node].keys(): # check if a path exists between starting and final node
            queue = []

        while queue:
            path_so_far = queue.popleft()  # Pop from the front of the queue
            current_node = path_so_far[-1]

            # If we reach the final node, append the path to results
            if current_node == final_node:
                paths.append(path_so_far)
                continue

            for neighbor in nx.all_neighbors(self.get_graph(), current_node):
                # Only explore neighbors that can potentially lead to the shortest path
                if  self.get_shortest_path_lengths()[neighbor][final_node] <  self.get_shortest_path_lengths()[current_node][final_node] and neighbor not in path_so_far:
                    queue.append(path_so_far + [neighbor])
        return paths

    def find_prime_mid_nodes(self, src_node:int) -> dict:
        '''
        Find all prime mid node in a graph from a source node.
        A prime mid node is a node that possesses at least two shortest paths of equal length to the source node.
        Combining two shortest paths should lead to a primitive ring.
        For each even rings, there is one prime-mid node
        While for each odd rings there are two linked prime mid nodes
        --------------------------------------------------------------------------------------------------------------
        Inputs:

            graph: networkx Graph object
            src_node: source node index
            shortest_dist_dict: dict[int] containing shortest paths length between the source node and each node in the graph
            max_ring_size: size of the largest ring to investigate

        Outputs:
            A dict with future primitive ring as key and prime mid node as values (for odd rings, it is a list of the linked prime mid nodes)
        
        '''
        # start_time = time.time()
        mid_nodes =  {x:[] for x in range(2, self.get_max_size() + 2)}
        src_shortest_paths = self.get_shortest_path_lengths()[src_node]
        for node, length in src_shortest_paths.items():
            if length <= self.get_max_size()//2 and length > 1:
                neighb_dist_list = [src_shortest_paths[neighbors] for neighbors in self.get_graph().neighbors(node)]
                if neighb_dist_list.count(length - 1) > 1:
                    mid_nodes[2*length].append(node)
                neighb_dist_dict = {neighbors:src_shortest_paths[neighbors] for neighbors in self.get_graph().neighbors(node)}
                for n, path in neighb_dist_dict.items():
                    if path == length and sorted([n, node]) not in mid_nodes[2*length + 1]:
                        mid_nodes[2*length + 1].append(sorted([n, node]))
        return mid_nodes
    
    def form_rings(self, src_node:int, prime_mid_node:dict) -> dict:
        '''
        Form a ring from prime mid nodes by combining two shortest paths.
        Exclude cases where two shortest paths share a node (other than initial or final)
        ------------------------------------------------------------------------------------
        Inputs:

            graph: networkx Graph object
            src_node: source node index
            shortest_dist_dict: dict[int] containing shortest paths length between the source node and each node in the graph
            prime_mid_node: dict[list] with ring size as keys and list of prime mid node as values
            max_ring_size: size of the largest ring to investigate

        Outputs:
            A dict with ring size keys and the shortest path to combine to form rings
        
        '''
        rings = {x:[] for x in range(2, self.get_max_size() + 2)}
        for ring_size, list_nodes in prime_mid_node.items():
            if ring_size%2 == 0:
                for nodes in list_nodes:
                    paths = self.find_all_shortest_paths(starting_node=src_node, final_node=nodes)
                    for poss_ring in combinations(paths, 2):
                        ring_a, ring_b = poss_ring
                        if len(set(ring_a[1:-1]) & set(ring_b[1:-1])) == 0:
                            rings[ring_size].append([ring_a, ring_b])
            elif ring_size%2 == 1:
                for lnked_nodes in list_nodes:
                    rings_a = self.find_all_shortest_paths(starting_node=src_node, final_node=lnked_nodes[0])
                    rings_b = self.find_all_shortest_paths(starting_node=src_node, final_node=lnked_nodes[1])
                    for ring_a in rings_a:
                        for ring_b in rings_b:
                            if len(set(ring_a[1:]) & set(ring_b[1:])) == 0:
                                rings[ring_size].append([ring_a, ring_b])
        return rings
    
    def exclude_non_primitive_rings(self, future_rings:dict) -> dict:
        '''
        Form primitive rings by combining two shortest paths.
        Exclude cases where there exist a shortcut between the two shortest paths
        ------------------------------------------------------------------------------------
        Inputs:

            potential_rings: dict[list] with ring size as keys and the shortest paths to combine to form rings
            shortest_dist_dict: dict[int] containing shortest paths length between the source node and each node in the graph
            max_ring_size: size of the largest ring to investigate

        Outputs:
            A dict with ring size keys and primitive rings as values
        
        '''
        missed_rings = {x:[] for x in range(2, self.get_max_size() + 2)}
        prim_rings =  {x:[] for x in range(2, self.get_max_size() + 2)}
        shortest_lengths= self.get_shortest_path_lengths()
        for ring_size, list_rings in future_rings.items():
            path_length = ring_size//2
            if len(list_rings) > 0:
                for ring_a, ring_b in list_rings:
                    shortcut = 0
                    if ring_size%2 == 0:
                        full_ring = ring_a + ring_b[1:-1][::-1]
                        for i in range(1, path_length):
                            check_node_a = full_ring[i]
                            check_node_b = full_ring[i + path_length]
                            if shortest_lengths[check_node_a][check_node_b] != path_length:
                                shortcut += 1
                                break
                        if shortcut == 0:
                            prim_rings[ring_size].append(full_ring)
                        elif shortcut == 1:
                            missed_rings[ring_size].append(full_ring)
                    elif ring_size%2 == 1:
                        full_ring = ring_a + ring_b[1:][::-1]
                        for i in range(1, path_length + 1):
                            check_node_a = full_ring[i]
                            check_node_b1 = full_ring[i + path_length]
                            check_node_b2 = full_ring[(1 + path_length + i)%ring_size]
                            check_a = shortest_lengths[check_node_a][check_node_b1] - path_length
                            check_b = shortest_lengths[check_node_a][check_node_b2] - path_length
                            if check_a + check_b != 0:
                                shortcut += 1
                                break
                        if shortcut == 0:
                            prim_rings[ring_size].append(full_ring)
                        elif shortcut == 1:
                            missed_rings[ring_size].append(full_ring)

        return prim_rings
    
    def search_primitive_ring(self):
        '''
        Search primitive rings following implementaion from Yuan et al. (https://doi.org/10.1016/S0927-0256(01)00256-7) in three steps:
        1) Find potential prime mid nodes
        2) Form rings from shortest paths between source nodes and prime nodes
        3) Exclude rings that contain shortcut (non-primitive)
        ------------------------------------------------------------------------------------------------------
        Outputs:
                A dict containing ring size as keys and list[list] of index forming rings
        '''
        start = time.time()
        all_rings = {x:[] for x in self.get_graph()}
        ring_count = {x:0 for x in range(2, self.get_max_size() + 2)}
        for src_node in self.src_nodes:
            prime_nodes = self.find_prime_mid_nodes(src_node=src_node)
            potential_rings = self.form_rings(src_node=src_node, prime_mid_node=prime_nodes)
            primitive_rings = self.exclude_non_primitive_rings(future_rings=potential_rings)
            all_rings[src_node] = primitive_rings
            for size in primitive_rings:
                ring_count[size] += len(primitive_rings[size])
        self.prim_rings = ring_count
        end = time.time()
        print(f"Computing primitive rings on {len(self.src_nodes)} atoms took in total {end - start} seconds")
        return all_rings
    
    def plot(self, color="blue", label_loc=0.005, normalized:bool=False):
        '''
        Plot a bar chart of ring sizes and occurence.
        -------------------------------------------------------------------
        Inputs:
                color: str a matploytlib color
                label_loc float locate label of ring occurence
        '''
        ring_per_si = np.array(list(self.prim_rings.values()))
        if normalized is True:
            ring_per_si = np.array(list(self.prim_rings.values()))/len(self.src_nodes)

        def addlabels(x,y, col, label_height):
            for i in range(len(x)):
                if y[i] - 0.00 > 0.05:
                    plt.text(x[i], y[i] + label_height*max(y), f"{y[i]:.2f}", ha = 'center', size=7, c=col, weight="bold")

        plt.bar(self.prim_rings.keys(), ring_per_si, color=color)
        addlabels(x=list(self.prim_rings.keys()), y=ring_per_si, col=color, label_height=label_loc)
        plt.xlabel("Ring size")


    
    
        