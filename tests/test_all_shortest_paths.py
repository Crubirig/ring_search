import networkx as nx
from searching_ring.Primitive_ring_search import Primitive_ring_search

def test_single_shortest_path():
    G = nx.Graph()
    G.add_edges_from([(1, 2), (2, 3)])
    sp_len = dict(nx.all_pairs_shortest_path_length(G))

    test_search = Primitive_ring_search(graph=G, max_size=4, distances=sp_len)

    result = test_search.find_all_shortest_paths(1, 3)
    assert result == [[1, 2, 3]]

def test_multiple_shortest_paths():
    G = nx.Graph()
    G.add_edges_from([(1, 2), (1, 4), (2, 3), (4, 3)])
    sp_len = dict(nx.all_pairs_shortest_path_length(G))

    test_search = Primitive_ring_search(graph=G, max_size=4, distances=sp_len)

    result = test_search.find_all_shortest_paths(1, 3)
    expected = [[1, 2, 3], [1, 4, 3]]
    assert sorted(result) == sorted(expected)

def test_no_path():
    G = nx.Graph()
    G.add_edges_from([(1, 2), (3, 4)])
    sp_len = dict(nx.all_pairs_shortest_path_length(G))
    test_search = Primitive_ring_search(graph=G, max_size=4, distances=sp_len)

    result = test_search.find_all_shortest_paths(1, 4)
    assert result == []

def test_start_equals_end():
    G = nx.Graph()
    G.add_edges_from([(1, 2), (2, 3)])
    sp_len = dict(nx.all_pairs_shortest_path_length(G))

    test_search = Primitive_ring_search(graph=G, max_size=4, distances=sp_len)

    result = test_search.find_all_shortest_paths(2, 2)
    assert result == [[2]]

