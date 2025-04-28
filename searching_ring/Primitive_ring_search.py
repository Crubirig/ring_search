import networkx as nx

from itertools import combinations
from collections import deque

class Primitive_ring_search():
    def __init__(self, graph:nx.Graph, max_size:int, distances:dict[dict]=None, src_nodes:list=[0]):
        self.graph = graph
        self.max_size = max_size
        self.src_nodes = src_nodes
        self.distances = distances if distances is not None else dict(
            nx.all_pairs_shortest_path_length(graph, cutoff=max_size//2 + 1))

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
        shortest_paths = self.get_shortest_path_lengths()[src_node]
        for node, length in shortest_paths.items():
            if length <= self.get_max_size()//2 and length > 1:
                neighb_dist_list = [shortest_paths[neighbors] for neighbors in self.get_graph().neighbors(node)]
                if neighb_dist_list.count(length - 1) > 1:
                    mid_nodes[2*length].append(node)
                neighb_dist_dict = {neighbors:shortest_paths[neighbors] for neighbors in self.get_graph().neighbors(node)}
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
                    for poss_ring in combinations(self.find_all_shortest_paths(starting_node=nodes, final_node=src_node), 2):
                        ring_a, ring_b = poss_ring
                        if len(set(ring_a[1:-1]) & set(ring_b[1:-1])) == 0:
                            rings[ring_size].append([ring_a, ring_b])
            elif ring_size%2 == 1:
                for lnked_nodes in list_nodes:
                    rings_a = self.find_all_shortest_paths(starting_node=lnked_nodes[0], final_node=src_node)
                    rings_b = self.find_all_shortest_paths(starting_node=lnked_nodes[1], final_node=src_node)
                    for ring_a in rings_a:
                        for ring_b in rings_b:
                            if len(set(ring_a[:-1]) & set(ring_b[:-1])) == 0:
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
        prim_rings =  {x:[] for x in range(2, self.get_max_size() + 2)}
        for ring_size, list_rings in future_rings.items():
            for rings in list_rings:
                ring_a, ring_b = rings
                shortcut = 0
                for i in range(1, len(ring_a) - 1):
                    if self.get_shortest_path_lengths()[ring_a[i]][ring_b[::-1][i]] != ring_size//2:
                        shortcut += 1
                        break
                if shortcut == 0:
                    if ring_size%2 == 0:
                        prim_rings[ring_size].append(ring_a[1:-1]+ring_b[::-1])
                    elif ring_size%2 == 1:
                        prim_rings[ring_size].append(ring_a[:-1]+ring_b[::-1])
        return prim_rings
    
    def single_node_search(self):
        all_rings = {x:[] for x in self.get_graph()}
        ring_count = {x:0 for x in range(2, self.get_max_size() + 2)}
        for src_node in self.src_nodes:
            prime_nodes = self.find_prime_mid_nodes(src_node=src_node)
            potential_rings = self.form_rings(src_node=src_node, prime_mid_node=prime_nodes)
            primitive_rings = self.exclude_non_primitive_rings(future_rings=potential_rings)
            all_rings[src_node] = primitive_rings
            for size in primitive_rings:
                ring_count[size] += len(primitive_rings[size])/len(self.src_nodes)
        return all_rings, ring_count
    
    


    
    
        