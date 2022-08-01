#!/bin/bash

#SBATCH --time=6:00:00
#SBATCH --mem=500M
#SBATCH --partition=batch-skl
#SBATCH --cpus-per-task=1
#SBATCH --ntasks=1
#SBATCH --output=edge-scale.out
make clean
make test16
make digraph16

SEED=123
VERTICES=40

for DEG in 1 2 3 4 6 8 11 14 18 23 30 40
do
    echo "degree: $DEG"
    ./digraph16 -qtf graphs/edge_scalability/cm${VERTICES}_${DEG}_0 -s $SEED

    EDGES=$(( $DEG * $VERTICES ))
    echo "edges: $EDGES"
    ./digraph16 -qtf graphs/edge_scalability/er${VERTICES}_${EDGES}_0 -s $SEED
done
