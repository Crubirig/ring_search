import networkx as nx
import time

from itertools import combinations

def find_all_shortest_paths(graph, starting_node:int, final_node:int, sortest_path_length:dict):
    #PREMIERE TENTATIVE
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

def primitive_ring_search(graph:nx.Graph, node_list:list=[0], max_ring_size:int=20):
    start_time = time.time()
    furthest_search = max_ring_size//2 + 1
    start_path = time.time()
    all_shortest_path_length = dict(nx.all_pairs_shortest_path_length(graph, cutoff=furthest_search))
    end_path = time.time()
    print(f"Computing all shortest path lengths took {end_path -start_path} seconds")
    primitive_ring_count = {x:0 for x in range(4, max_ring_size + 1)}
    for src_node in node_list:
        prime_mid_node = {x:[] for x in range(4, max_ring_size + 1)}
        future_rings  = {x:[] for x in range(4, max_ring_size + 1)}
        primitive_ring =  {x:[] for x in range(4, max_ring_size + 1)}

        shortest_paths = all_shortest_path_length[src_node]
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
        # end_time = time.time()
        # print(f"It took {end_time - start_time} seconds to find prime mid nodes")
        # start_time = time.time()
        for ring_size, list_nodes in prime_mid_node.items():
            if ring_size%2 == 0:
                for nodes in list_nodes:
                    for poss_ring in combinations(find_all_shortest_paths(graph, nodes, src_node, all_shortest_path_length), 2):
                        ring_a, ring_b = poss_ring
                        if len(set(ring_a[1:-1]) & set(ring_b[1:-1])) == 0:
                            future_rings[ring_size].append([ring_a, ring_b])
            elif ring_size%2 == 1:
                for lnked_nodes in list_nodes:
                    rings_a = find_all_shortest_paths(graph, lnked_nodes[0], src_node, all_shortest_path_length)
                    rings_b = find_all_shortest_paths(graph, lnked_nodes[1], src_node, all_shortest_path_length)
                    for ring_a in rings_a:
                        for ring_b in rings_b:
                            if len(set(ring_a[:-1]) & set(ring_b[:-1])) == 0:
                                future_rings[ring_size].append([ring_a, ring_b])
        # end_time = time.time()
        # print(f"It took {end_time - start_time} seconds to form rings")
        # start_time = time.time()
        for ring_size, list_rings in future_rings.items():
            for rings in list_rings:
                ring_a, ring_b = rings
                shortcut = 0
                for i in range(1, len(ring_a) - 1):
                    if all_shortest_path_length[ring_a[i]][ring_b[::-1][i]] !=  ring_size//2:
                        shortcut += 1
                        break
                if shortcut == 0:
                    if ring_size%2 == 0:
                        primitive_ring[ring_size].append(ring_a[1:-1]+ring_b[::-1])
                        primitive_ring_count[ring_size] += 1
                    elif ring_size%2 == 1:
                        primitive_ring[ring_size].append(ring_a[:-1]+ring_b[::-1])
                        primitive_ring_count[ring_size] += 1
        # end_time = time.time()
        # print(f"It took {end_time - start_time} seconds to exclude non-primitive rings")
    end_time = time.time()
    print(f"Computing all primitive rings up to {max_ring_size}-rings on {len(node_list)} atoms took {end_time - start_time} seconds")
    print(primitive_ring_count)