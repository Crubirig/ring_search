import networkx as nx
import time

from itertools import combinations

import networkx as nx
from collections import deque

def find_all_shortest_paths(graph: nx.Graph, starting_node: int, final_node: int, shortest_path_length: dict) -> list[list]:
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

    while queue:
        path_so_far = queue.popleft()  # Pop from the front of the queue
        current_node = path_so_far[-1]

        # If we reach the final node, append the path to results
        if current_node == final_node:
            paths.append(path_so_far)
            continue

        for neighbor in nx.all_neighbors(graph, current_node):
            # Only explore neighbors that can potentially lead to the shortest path
            if shortest_path_length[neighbor][final_node] < shortest_path_length[current_node][final_node] and neighbor not in path_so_far:
                queue.append(path_so_far + [neighbor])

    return paths

def old_find_all_shortest_paths(graph:nx.Graph, starting_node:int, final_node:int, sortest_path_length:dict) -> list[list]:
    '''
    Find all shortest paths between two nodes in a graph knowing the length of those paths. 
    This is much faster than when the length is not known.
    -------------------------------------------------------------------------------------------------------
    Inputs: 
    
            graph: nteworkX graph
            starting node: starting node index
            final node: target node index
            shortest_path_length: dict[dict[int]] containing shortest paths length for each node pairs

    Outputs:
            A list of every shortest paths between starting and target node

    '''
    queue = [[starting_node]]
    path = []
    while len(queue) > 0:
        for i in queue:
            nodcrt = i
            for j in nx.all_neighbors(graph, nodcrt[-1]):
                if sortest_path_length[j][final_node] < sortest_path_length[nodcrt[-1]][final_node] and j not in nodcrt:
                    queue.append(nodcrt + [j])
                    if sortest_path_length[j][final_node] == 0:
                        path.append(nodcrt + [j])
            queue.remove(nodcrt)
    return path

def find_prime_mid_nodes(src_node:int, shortest_dist_dict:dict[dict[list]], max_ring_size:int, graph:nx.Graph):
    shortest_paths = shortest_dist_dict[src_node]
    prime_mid_node = {x:[] for x in range(2, max_ring_size + 1)}
    # start_time = time.time()
    for node, length in shortest_paths.items():
        if length <= max_ring_size//2 and length > 1:
            neighb_dist_list = [shortest_paths[neighbors] for neighbors in graph.neighbors(node)]
            if neighb_dist_list.count(length - 1) > 1:
                prime_mid_node[2*length].append(node)
            neighb_dist_dict = {neighbors:shortest_paths[neighbors] for neighbors in graph.neighbors(node)}
            for n, path in neighb_dist_dict.items():
                if path == length and sorted([n, node]) not in prime_mid_node[2*length + 1]:
                    prime_mid_node[2*length + 1].append(sorted([n, node]))
    return prime_mid_node

def form_rings(prime_mid_node:dict, src_node:int, shortest_dist_dict:dict[dict[list]], max_ring_size:int, graph:nx.Graph):
    future_rings  = {x:[] for x in range(2, max_ring_size + 1)}
    for ring_size, list_nodes in prime_mid_node.items():
        if ring_size%2 == 0:
            for nodes in list_nodes:
                for poss_ring in combinations(find_all_shortest_paths(graph, nodes, src_node, shortest_dist_dict), 2):
                    ring_a, ring_b = poss_ring
                    if len(set(ring_a[1:-1]) & set(ring_b[1:-1])) == 0:
                        future_rings[ring_size].append([ring_a, ring_b])
        elif ring_size%2 == 1:
            for lnked_nodes in list_nodes:
                rings_a = find_all_shortest_paths(graph, lnked_nodes[0], src_node, shortest_dist_dict)
                rings_b = find_all_shortest_paths(graph, lnked_nodes[1], src_node, shortest_dist_dict)
                for ring_a in rings_a:
                    for ring_b in rings_b:
                        if len(set(ring_a[:-1]) & set(ring_b[:-1])) == 0:
                            future_rings[ring_size].append([ring_a, ring_b])
    return future_rings

def exclude_non_primitive_rings(potential_rings:dict, max_ring_size:int, shortest_dist_dict:dict[dict[list]]):
    primitive_ring = {x:[] for x in range(2, max_ring_size + 1)}
    for ring_size, list_rings in potential_rings.items():
        for rings in list_rings:
            ring_a, ring_b = rings
            shortcut = 0
            for i in range(1, len(ring_a) - 1):
                if shortest_dist_dict[ring_a[i]][ring_b[::-1][i]] !=  ring_size//2:
                    shortcut += 1
                    break
            if shortcut == 0:
                if ring_size%2 == 0:
                    primitive_ring[ring_size].append(ring_a[1:-1]+ring_b[::-1])
                    # primitive_ring_count[ring_size] += 1
                elif ring_size%2 == 1:
                    primitive_ring[ring_size].append(ring_a[:-1]+ring_b[::-1])
                    # primitive_ring_count[ring_size] += 1
    return primitive_ring

def find_primitive_ring(all_path_length:dict[dict[list]], graph:nx.Graph, src_node:int, max_dist:int):
    nodes = find_prime_mid_nodes(src_node=src_node, shortest_dist_dict=all_path_length, max_ring_size=max_dist, graph=graph)

    rings = form_rings(prime_mid_node=nodes,src_node=src_node,shortest_dist_dict=all_path_length,max_ring_size=max_dist, graph=graph)

    primitive = exclude_non_primitive_rings(potential_rings=rings, max_ring_size=max_dist, shortest_dist_dict=all_path_length)
    return primitive

def primitive_ring_search(graph:nx.Graph, node_list:list=None, max_ring_size:int=20, verbose=True) -> dict:
    if node_list is None:
        node_list = graph.nodes
        print("Perform ring search on the entire graph.")
    start_time = time.time()
    furthest_search = max_ring_size//2 + 1
    start_path = time.time()
    all_paths = dict(nx.all_pairs_shortest_path_length(graph, cutoff=furthest_search))
    end_path = time.time()
    print(f"Computing all shortest path lengths took {end_path -start_path} seconds")
    all_rings = {x:dict() for x in node_list}
    for src_node in node_list:
        all_rings[src_node] = find_primitive_ring(src_node=src_node, all_path_length=all_paths, graph=graph, max_dist=max_ring_size)
    end_time = time.time()
    print(f"Computing all primitive rings up to {max_ring_size}-rings on {len(node_list)} atoms took {end_time - start_time} seconds")
    # print(primitive_ring_count)
    if verbose is True:
        primitive_ring_count = {x:0 for x in range(2, max_ring_size + 1)}
        for src_node in node_list:
            for size, rings in all_rings[src_node].items():
                primitive_ring_count[size] += len(rings)
        print(primitive_ring_count)
    return all_rings