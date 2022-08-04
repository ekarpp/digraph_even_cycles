#!/bin/bash

#SBATCH --time=??:00:00
#SBATCH --mem=1G
#SBATCH --partition=batch-skl
#SBATCH --cpus-per-task=40
#SBATCH --ntasks=1
#SBATCH --output=results/topo/algo-topology-par.out

hostname
uname -a
cat /etc/*release
g++ --version

make clean
make digraph-tests16
make digraph16

make test16

SEED=123

do_topology() {
    for vertices in 16 24 32 40 48 56
    do
        echo "vertices: ${vertices}"
        for i in {0..5}
        do
            graph=graphs/topology_variance/${1}${vertices}_*_${i}
            ./digraph16-PAR -qtf $graph
        done
    done
}

do_cycle_or_complete() {
    for vertices in 16 24 32 40 48 56
    do
        echo "vertices: ${vertices}"
        for i in {0..5}
        do
            graph=graphs/topology_variance/${1}${vertices}
            ./digraph16-PAR -qtf $graph
        done
}

echo "ERDOS"
do_topology "er"

echo "CONFIG MODEL"
do_topology "cm"

echo "CYCLE"
do_cycle_or_complete c

echo "COMPLETE"
do_cycle_or_complete k
