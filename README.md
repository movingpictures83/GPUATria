# GPUATria
# Language: CUDA
# Input: CSV (network)
# Output: NOA (central nodes and centrality values)

PluMA plugin to run the Ablatio Triadum (ATria, Cickovski et al 2015, 2017) algorithm on the GPU.
ATria can find the most important or central nodes in a signed and weighted network.
This plugin accepts input in the form of a CSV file with rows and columns corresponding to nodes,
and entry (i, j) the weight from node i to node j.

The plugin will then produce an output NOde Attribute file containing a simple table
that maps node names to centrality values.  This file is then readable by Cytoscape,
with Centrality becoming an attribute of every node in the network (in addition to its Rank).
These attributes can then be used for further downstream analysis or visualization.
