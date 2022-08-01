#!/bin/bash

#SBATCH --time=3:00:00
#SBATCH --mem=1G
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

for E in 1 2 3 4 5 7 9 11 14 17 21 26 32 39
do
    echo "edges: $E"
    for i in {0..5}
    do
        echo cm40_${E}_${i}
        ./digraph16 -qtf graphs/edge_scalability/cm40_${E}_${i} -s $SEED
    done
done
