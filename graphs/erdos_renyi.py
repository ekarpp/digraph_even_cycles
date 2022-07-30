#!/usr/bin/env python3
import random
import itertools
random.seed(123)

# random graph of V vertices and E edges
# edges are chosen uniformly at random
# from all possibilities
def erdos_renyi(V, E, i):
    vertices = [x for x in range(V)]
    pairs = list(itertools.combinations(vertices, 2))
    random.shuffle(pairs)
    edges = pairs[:E]
    adj = [[] for _ in range(V)]
    for u,v in edges:
        adj[u].append(str(v))
        adj[v].append(str(u))
    grph = "\n".join([" ".join(l) for l in adj]) + "\n"
    with open(f"erdos_renyi/er{V}_{E}_{i}", "w") as f:
        f.write(grph)

import sys

V = int(sys.argv[1])
E = int(sys.argv[2])
n = int(sys.argv[3])

for i in range(n):
    erdos_renyi(V, E, i)
