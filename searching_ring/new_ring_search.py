import networkx as nx
import time
from itertools import combinations

def primitive_ring_search(graph:nx.Graph, node_list:list=[0], max_ring_size:int=20):
    furthest_search = max_ring_size//2 + 1
    start_time = time.time()
    all_shortest_path_length = dict(nx.all_pairs_shortest_path_length(graph, cutoff=furthest_search))
    end_time = time.time()
    print(f"Computing all shortest paths took {end_time -start_time} seconds")
    primitive_ring_count = {x:0 for x in range(4, max_ring_size + 1)}
    for src_node in node_list:
        prime_mid_node = {x:[] for x in range(4, max_ring_size + 1)}
        future_rings  = {x:[] for x in range(4, max_ring_size + 1)}
        primitive_ring =  {x:[] for x in range(4, max_ring_size + 1)}

        shortest_paths = all_shortest_path_length[src_node]
        start_time = time.time()
        for node, length in shortest_paths.items():
            if length <= max_ring_size//2 and length > 1:
                neighb_dist_list = [shortest_paths[neighbors] for neighbors in graph.neighbors(node)]
                if neighb_dist_list.count(length - 1) > 1:
                    prime_mid_node[2*length].append(node)
                neighb_dist_dict = {neighbors:shortest_paths[neighbors] for neighbors in graph.neighbors(node)}
                for n, path in neighb_dist_dict.items():
                    if path == length and sorted([n, node]) not in prime_mid_node[2*length + 1]:
                        prime_mid_node[2*length + 1].append(sorted([n, node]))
        end_time = time.time()
        print(f"It took {end_time - start_time} seconds to find prime mid nodes")
        start_time = time.time()
        for ring_size, list_nodes in prime_mid_node.items():
            if ring_size%2 == 0:
                for nodes in list_nodes:
                    for poss_ring in combinations(nx.all_shortest_paths(graph, nodes, src_node), 2):
                        ring_a, ring_b = poss_ring
                        if len(set(ring_a[1:-1]) & set(ring_b[1:-1])) == 0:
                            future_rings[ring_size].append([ring_a, ring_b])
            elif ring_size%2 == 1:
                for lnked_nodes in list_nodes:
                    rings_a = nx.all_shortest_paths(graph, lnked_nodes[0], src_node)
                    rings_b = nx.all_shortest_paths(graph, lnked_nodes[1], src_node)
                    for ring_a in rings_a:
                        for ring_b in rings_b:
                            if len(set(ring_a[:-1]) & set(ring_b[:-1])) == 0:
                                future_rings[ring_size].append([ring_a, ring_b])
        end_time = time.time()
        print(f"It took {end_time - start_time} seconds to form rings")
        start_time = time.time()
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
        end_time = time.time()
        print(f"It took {end_time - start_time} seconds to exclude non-primitive rings")
        print(src_node)
                        
    print(primitive_ring_count)