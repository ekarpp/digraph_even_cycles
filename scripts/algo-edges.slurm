#!/bin/bash

#SBATCH --time=6:00:00
#SBATCH --mem=500M
#SBATCH --partition=batch-skl
#SBATCH --cpus-per-task=1
#SBATCH --ntasks=1
#SBATCH --output=edge-scale.out
module load gcc/11.2.0

hostname
uname -a
cat /etc/*release
g++ --version

make clean
make test16
make digraph16

SEED=123
VERTICES=40

for DEG in 1 2 3 4 6 8 11 14 18 23 30 40
do
    echo "degree: $DEG"
    for i in {0..5}
    do
        ./digraph16 -qtf graphs/edge_scalability/cm${VERTICES}_${DEG}_${i} -s $SEED
    done

    EDGES=$(( $DEG * $VERTICES ))
    echo "edges: $EDGES"
    for i in {0..5}
    do
        ./digraph16 -qtf graphs/edge_scalability/er${VERTICES}_${EDGES}_${i} -s $SEED
    done
done
