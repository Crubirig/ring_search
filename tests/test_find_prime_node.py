import networkx as nx

from searching_ring.new_ring_search import find_prime_mid_nodes, form_rings, exclude_non_primitive_rings

def test_find_prime_mid_nodes():
    edges = [(i, (i+1+6)%6) for i in range(6)]
    five_cycle = [(4, 6), (6, 7), (7, 8), (8, 5)]
    edges += five_cycle
    G = nx.Graph(edges)

    expected_result = {2: [], 3: [], 4: [], 5: [], 6: [4], 7: [], 8: [], 9: [[6, 7]], 10: []}

    shortest_dist_dict = dict(nx.all_pairs_shortest_path_length(G))
    result = find_prime_mid_nodes(graph=G, src_node=1, shortest_dist_dict=shortest_dist_dict, max_ring_size=9)

    assert result == expected_result

def test_simple_form_rings():
    edges = [(i, (i+1+6)%6) for i in range(6)]
    five_cycle = [(4, 6), (6, 7), (7, 8), (8, 5)]
    edges += five_cycle
    G = nx.Graph(edges)

    shortest_dist_dict = dict(nx.all_pairs_shortest_path_length(G))
    mid_node_list = {2: [], 3: [], 4: [], 5: [], 6: [4], 7: [], 8: [], 9: [[6, 7]], 10: []}

    result = form_rings(graph=G, src_node=1, prime_mid_node=mid_node_list, shortest_dist_dict=shortest_dist_dict, max_ring_size=9)

    expected_result = {2: [], 3: [], 4: [], 5: [], 6: [[[4, 3, 2, 1], [4, 5, 0, 1]]], 7: [], 8: [], 9: [[[6, 4, 3, 2, 1], [7, 8, 5, 0, 1]]], 10: []}

    assert result == expected_result

def test_form_rings_exclude_shared_node():
    edges = [(i, (i+1+6)%6) for i in range(6)]
    five_cycle = [(4, 6), (6, 7), (7, 8), (8, 5)]
    edges += five_cycle
    G = nx.Graph(edges)

    shortest_dist_dict = dict(nx.all_pairs_shortest_path_length(G))
    mid_node_list = {2: [], 3: [], 4: [], 5: [], 6: [0], 7: [[7, 8]], 8: [], 9: [], 10: []}

    result = form_rings(graph=G, src_node=3, prime_mid_node=mid_node_list, shortest_dist_dict=shortest_dist_dict, max_ring_size=9)
    
    expected_result = {2: [], 3: [], 4: [], 5: [], 6: [[[0, 1, 2, 3], [0, 5, 4, 3]]], 7: [], 8: [], 9: [], 10: []}

    assert result == expected_result

def test_exclude_primitive_rings_nothing_to_exclude():
    edges = [(i, (i+1+6)%6) for i in range(6)]
    five_cycle = [(4, 6), (6, 7), (7, 8), (8, 5)]
    edges += five_cycle
    G = nx.Graph(edges)

    shortest_dist_dict = dict(nx.all_pairs_shortest_path_length(G))
    future_rings = {2: [], 3: [], 4: [], 5: [], 6: [[[0, 1, 2, 3], [0, 5, 4, 3]]], 7: [], 8: [], 9: [], 10: []}

    result = exclude_non_primitive_rings(potential_rings=future_rings, max_ring_size=9, shortest_dist_dict=shortest_dist_dict)
    
    expected_result = {2: [], 3: [], 4: [], 5: [], 6: [[1, 2, 3, 4, 5, 0]], 7: [], 8: [], 9: [], 10: []}

    assert result == expected_result

def test_exclude_primitive_rings_shortcut():
    edges = [(i, (i+1+6)%6) for i in range(6)]
    five_cycle = [(4, 6), (6, 7), (7, 8), (8, 5)]
    edges += five_cycle
    G = nx.Graph(edges)

    shortest_dist_dict = dict(nx.all_pairs_shortest_path_length(G))
    future_rings = {2: [], 3: [], 4: [], 5: [], 6: [[[4, 3, 2, 1], [4, 5, 0, 1]]], 7: [], 8: [], 9: [[[6, 4, 3, 2, 1], [7, 8, 5, 0, 1]]], 10: []}

    result = exclude_non_primitive_rings(potential_rings=future_rings, max_ring_size=9, shortest_dist_dict=shortest_dist_dict)
    
    expected_result = {2: [], 3: [], 4: [], 5: [], 6: [[3, 2, 1, 0, 5, 4]], 7: [], 8: [], 9: [], 10: []}

    assert result == expected_result