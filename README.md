# ring_search
Python script to compute primitive rings of silica structure

## Overview
This project provides tools to analyze silica structures and compute primitive rings using graph theory. It leverages the `networkx` library for graph operations and the `ase` library for handling atomic structures.

## Key Modules

### `Primitive_ring_search.py`
This module defines the `Primitive_ring_search` class, which implements a graph-based algorithm to find primitive rings in a molecular structure. The algorithm is based on the work of Yuan and Cormack.

#### Features:
- **Shortest Path Calculation**: Efficiently computes all shortest paths between nodes in the graph.
- **Primitive Ring Detection**: Identifies primitive rings by combining shortest paths and excluding shortcuts.
- **Visualization**: Provides a method to plot the distribution of ring sizes.

#### Example Usage:
```python
from searching_ring.Primitive_ring_search import Primitive_ring_search
import networkx as nx

# Create a graph
graph = nx.Graph()
graph.add_edges_from([(0, 1), (1, 2), (2, 3), (3, 4), (4, 5), (5, 0)])

# Initialize the search
search = Primitive_ring_search(graph=graph, max_size=6)

# Perform the search
primitive_rings = search.search_primitive_ring()

# Plot the results
search.plot()