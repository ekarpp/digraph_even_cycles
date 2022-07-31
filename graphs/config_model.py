#!/usr/bin/env python3
import random

def create_adja(heads, tails, n):
    adj = [[] for _ in range(n)]
    while heads:
        idx = random.randrange(0, len(heads))
        h = heads.pop(idx)
        idx = random.randrange(0, len(tails))
        t = tails.pop(idx)
        adj[t].append(h)
    return adj

def adja_ok(adj):
    for i in range(len(adj)):
        if i in adj[i] or len(adj[i]) != len(set(adj[i])):
            return False
    return True

def adja_ok_dups(adj):
    for i in range(len(adj)):
        if i in adj[i]:
            return False
    return True

def cfg_model(n, d, m, path, dups):
    W = []
    for i in range(n):
        for j in range(d):
            W.append(i)
    G = create_adja(W.copy(), W.copy(), n)
    test_adj = adja_ok
    if dups:
        test_adj = adja_ok_dups
    while not test_adj(G):
        G = create_adja(W.copy(), W.copy(), n)
    graph = "\n".join(
        [" ".join([str(x) for x in l]) for l in G]
    ) + "\n"
    with open(f"{path}/cm{n}_{d}_{m}", "w") as f:
        f.write(graph)

import sys

n = int(sys.argv[1])
d = int(sys.argv[2])
m = int(sys.argv[3])
p = sys.argv[4]
dups = len(sys.argv) > 5

if n*d % 2 == 0:
    for i in range(m):
        cfg_model(n, d, i, p, dups)
