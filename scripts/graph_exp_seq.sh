#!/bin/bash

#SBATCH --time=6:00:00
#SBATCH --mem=1G
#SBATCH --partition=batch-hsw
#SBATCH --cpus-per-task=24
#SBATCH --ntasks=1
#SBATCH --output=results/graph_exp.out

hostname
uname -a
cat /etc/*release
g++ --version

make clean
make objects -j24
make digraph
make test

SEED=123

edge_scale_config() {
    VERTICES=40
    for DEG in 1 2 3 4 6 8 11 14 18 23 30 39
    do
        echo "config degree: $DEG"
        for i in {0..5}
        do
            graph=graphs/${1}/cm${vertices}_${DEG}_${i}
            ${2} -qtf ${graph} -s $SEED
        done
    done
}

edge_scale_er() {
    VERTICES=40
    for DEG in 1 2 3 4 6 8 11 14 18 23 30 39
    do
        echo "erdos degree: $DEG"
        EDGES=$(( $DEG * $VERTICES ))
        for i in {0..5}
        do
            graph=graphs/${1}/er${vertices}_${EDGES}_${i}
            ${2} -qtf ${graph} -s $SEED
        done
    done
}

do_topology() {
    for VERTICES in 16 24 32 40 48
    do
        echo "${2} VERTICES: ${VERTICES}"
        for i in {0..5}
        do
            graph=graphs/${1}/${2}${vertices}_*_${i}
            ${3} -qtf $graph -s $SEED
        done
    done
}

do_cycle_or_complete() {
    for VERTICES in 16 24 32 40 48
    do
        echo "${2} VERTICES: ${VERTICES}"
        for i in {0..5}
        do
            graph=graphs/${1}/${2}${vertices}
            ${3} -qtf $graph -s $SEED
        done
    done
}

# TOPOLOGY
BIN="./digraph16"
FOLDER="topology_variance"
RES_PTH="results/topo/seq"
mkdir -p ${RES_PTH}
do_topology $FOLDER "er" $BIN > ${RES_PTH}/er &
do_topology $FOLDER "cm" $BIN > ${RES_PTH}/cm &
do_cycle_or_copmlete $FOLDER "c" $BIN > ${RES_PTH}/c &
do_cycle_or_copmlete $FOLDER "k" $BIN > ${RES_PTH}/k &

# EDGE SCALABILITY
BIN="./digraph16"
FOLDER="edge_scalability"
RES_PTH="results/edge"
mkdir -p ${RES_PTH}
edge_scale_config $FOLDER $BIN > ${RES_PTH}/cm &
edge_scale_er $FOLDER $BIN > ${RES_PTH}/er &

# FF COMPARE
FOLDER="ff_comparison"
RES_PTH="results/ff_compare"
mkdir -p ${RES_PTH}
for b in 0 16 32
do
    do_topology $FOLDER "cm" "./digraph${b}" > ${RES_PTH}/digraph${b} &
done

# VECTOR PERF
FOLDER="vector_performance"
RES_PTH="results/vecotr_perf"
mkdir -p ${RES_PTH}
do_topology $FOLDER "cm" "./digraph16" > ${RES_PTH}/vec &
do_topology $FOLDER "cm" "./digraph16-NOVEC" > ${RES_PTH}/no_vec &

wait
zip -r $(date '+%Y-%m-%d')_results.zip results
