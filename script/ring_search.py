import networkx as nx
import time
import numpy as np

def ring_counting(T_atoms:list, silica_graph, max_ring=14):
    Total_ring_amount, primitive_ring = rings_search(SI_atoms=T_atoms, graph=silica_graph, Max_ring_size=max_ring)
    for i in np.arange(2, max_ring):
        Total_ring_amount[i] = 0
        for SI_index in T_atoms:                    
            Total_ring_amount[i] += len(primitive_ring[f"SI{SI_index}"][i])
    ring_per_si = dict()
    for i in Total_ring_amount:
        print(f"{Total_ring_amount[i]/len(T_atoms):.1f} {i}-membered ring(s) per Si atom")
        ring_per_si[i] = Total_ring_amount[i]/len(T_atoms)
    return primitive_ring, ring_per_si
    

def find_shortest_path(graph, starting_node:int, final_node:int, sortest_path_length:dict):
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

def rings_search(SI_atoms:list, graph, Max_ring_size:int=14):
    #Ring search algorithm implementaion based on Yuan and Cormack https://doi.org/10.1016/S0927-0256(01)00256-7
    #Algorithm using networkx library specialized in graph theory
    # Can be decomposed in three stpes
    print("Starting graph building")
    time_to_graph = time.time()
    all_shortest_path_length = dict(nx.all_pairs_shortest_path_length(graph, cutoff=int(Max_ring_size/2) + 2))
    graph_finished = time.time()
    print(f"It took {graph_finished - time_to_graph} seconds to compute shortest_path_length")
    shortest_path = {}
    primitive_ring = {}
    Total_ring_amount = {}
    for SI_index in SI_atoms:
        primitive_ring[f"SI{SI_index}"] = {}
        shortest_path[SI_index] = {}
        # STEP 1:Searching for potential prime_mid nodes among the lvldist (list of all shortest paths)
        # A prime mid node is the furthest node from the source node ("the first node of the ring") in a primitive ring
        # There is one prime mid node in a 2n ring and two in a 2n+1 ring
        prime_mid_node = {}
        for i in np.arange(1, int(Max_ring_size/2)):
            prime_mid_node[i*2] = []
            prime_mid_node[i*2 + 1] = []
            for k in all_shortest_path_length[SI_index]:
                if all_shortest_path_length[SI_index][k] == i:
                    just_below = 0
                    similar = 0
                    SI_neighbour = nx.all_neighbors(graph, k)
                    linked_nodes = [k]
                    for l in SI_neighbour:
                        if all_shortest_path_length[SI_index][l] == all_shortest_path_length[SI_index][k] - 1:
                            just_below +=1
                        elif all_shortest_path_length[SI_index][l] == all_shortest_path_length[SI_index][k]:
                            linked_nodes.append(l)
                            similar +=1
                    if just_below > 1:
                        prime_mid_node[i*2].append(k)
                        shortest_path[SI_index][k] = find_shortest_path(graph, SI_index, k, all_shortest_path_length)
                    if similar != 0:
                        shortest_path[SI_index][k] = find_shortest_path(graph, SI_index, k, all_shortest_path_length)
                        if sorted(linked_nodes) not in prime_mid_node[i*2 + 1]:
                            prime_mid_node[i*2 + 1].append(sorted(linked_nodes)) 
        #STEP 2: Combine two prime_mid_nodes to form a potential primitive ring
        short_path = {}
        for i in np.arange(2, Max_ring_size):
            short_path[i] = []
            if i%2 == 1:
                for j in prime_mid_node[i]:
                    half_path_1 = shortest_path[SI_index][j[0]]
                    half_path_2 = shortest_path[SI_index][j[1]]
                    for k in half_path_1:
                        for l in half_path_2:
                            if not any(a in l[1:-1] for a in k[1:-1]):
                                short_path[i].append(k + l[::-1])
            elif i%2 ==0:
                for j in prime_mid_node[i]:
                    for k in range(len(shortest_path[SI_index][j])):
                        for l in range(k, len(shortest_path[SI_index][j])):
                            if not any (a in shortest_path[SI_index][j][l][1:-1] for a in shortest_path[SI_index][j][k][1:-1]):
                                short_path[i].append(shortest_path[SI_index][j][k][:-1] + shortest_path[SI_index][j][l][::-1])
        #STEP 3:Among all primitive ring, delete those who contain shortcut
        for i in range(1, int(Max_ring_size/2)):
            primitive_ring[f"SI{SI_index}"][i*2] = []
            primitive_ring[f"SI{SI_index}"][i*2+1] = []
            for j in range(len(short_path[i*2])):
                short_cut = 0
                for k in range(1, i):
                    if all_shortest_path_length[short_path[i*2][j][k]][short_path[i*2][j][k+i]] < i:
                        short_cut += 1
                        break
                if short_cut == 0:
                    primitive_ring[f"SI{SI_index}"][i*2].append(short_path[i*2][j])
            for j in range(len(short_path[i*2 + 1])):
                short_cut = 0
                for k in range(1, i):
                    if all_shortest_path_length[short_path[i*2 + 1][j][k]][short_path[i*2 + 1][j][k+i]] < i or all_shortest_path_length[short_path[i*2 + 1][j][k]][short_path[i*2 + 1][j][k+i + 1]] < i:
                        short_cut += 1
                        break
                if short_cut == 0:
                    primitive_ring[f"SI{SI_index}"][i*2 + 1].append(short_path[i*2 + 1][j])
    return([Total_ring_amount, primitive_ring])


