#!/usr/bin/env python3
import random
import itertools

# random graph of V vertices and E edges
# edges are chosen uniformly at random
# from all possibilities
def erdos_renyi(V, E):
    vertices = [x for x in range(V)]
    pairs = list(itertools.combinations(vertices, 2))
    random.shuffle(pairs)
    edges = pairs[:E]
    adj = [[] for _ in range(V)]
    for u,v in edges:
        adj[u].append(str(v))
        adj[v].append(str(u))
    grph = "\n".join([" ".join(l) for l in adj]) + "\n"
    with open(f"er_{V}_{E}", "w") as f:
        f.write(grph)

import sys

V = int(sys.argv[1])
E = int(sys.argv[2])

erdos_renyi(V, E)
