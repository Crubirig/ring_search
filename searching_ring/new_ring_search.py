import networkx as nx
from itertools import combinations

def primitive_ring_search(graph:nx.Graph, src_node:int=0, max_ring_size:int=10):
    src_node = 0
    max_ring_size = 10
    furthest_search = max_ring_size//2 + 1
    all_shortest_path_length = dict(nx.all_pairs_shortest_path_length(graph, cutoff=furthest_search))
    primitive_ring_count = {x:0 for x in range(4, max_ring_size)}

    prime_mid_node = {x:[] for x in range(4, max_ring_size)}
    future_rings  = {x:[] for x in range(4, max_ring_size)}
    primitive_ring =  {x:[] for x in range(4, max_ring_size)}

    shortest_paths = all_shortest_path_length[src_node]
    for ring_size in range(4, max_ring_size):
        for node, length in shortest_paths.items():
            if length == ring_size//2:
                if ring_size%2 == 0:
                    neighb_dist = [all_shortest_path_length[neighbors][src_node] for neighbors in graph.neighbors(node)]
                    if neighb_dist.count(ring_size//2 - 1) > 1:
                        prime_mid_node[ring_size].append(node)
                else:
                    neighb_dist = {neighbors:all_shortest_path_length[neighbors][src_node] for neighbors in graph.neighbors(node)}
                    for n, path in neighb_dist.items():
                        if path == length and sorted([n, node]) not in prime_mid_node[ring_size]:
                            prime_mid_node[ring_size].append(sorted([n, node]))

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
    print(primitive_ring)
    print(primitive_ring_count)