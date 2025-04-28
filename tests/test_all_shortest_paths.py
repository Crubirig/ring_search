import networkx as nx
from searching_ring.new_ring_search import find_all_shortest_paths

def test_single_shortest_path():
    G = nx.Graph()
    G.add_edges_from([(1, 2), (2, 3)])
    sp_len = dict(nx.all_pairs_shortest_path_length(G))

    result = find_all_shortest_paths(G, 1, 3, sp_len)
    assert result == [[1, 2, 3]]

def test_multiple_shortest_paths():
    G = nx.Graph()
    G.add_edges_from([(1, 2), (1, 4), (2, 3), (4, 3)])
    sp_len = dict(nx.all_pairs_shortest_path_length(G))

    result = find_all_shortest_paths(G, 1, 3, sp_len)
    expected = [[1, 2, 3], [1, 4, 3]]
    assert sorted(result) == sorted(expected)

def test_no_path():
    G = nx.Graph()
    G.add_edges_from([(1, 2), (3, 4)])
    sp_len = dict(nx.all_pairs_shortest_path_length(G))

    result = find_all_shortest_paths(G, 1, 4, sp_len)
    assert result == []

def test_start_equals_end():
    G = nx.Graph()
    G.add_edges_from([(1, 2), (2, 3)])
    sp_len = dict(nx.all_pairs_shortest_path_length(G))

    result = find_all_shortest_paths(G, 2, 2, sp_len)
    assert result == [[2]]

