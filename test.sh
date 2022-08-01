#!/bin/bash

SEED=123

for E in 1 2 3 4 5 7 9 11 14 17 21 26 32 39
do
    echo "edges: $E"
    echo cm40_${E}_0
    ./digraph16 -qtf graphs/edge_scalability/cm40_${E}_0 -s $SEED
done
