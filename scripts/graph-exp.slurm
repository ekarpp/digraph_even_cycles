#!/bin/bash

# time not accurate
#SBATCH --time=4:00:00
#SBATCH --mem=1G
#SBATCH --partition=short-hsw
#SBATCH --cpus-per-task=24
#SBATCH --ntasks=1
#SBATCH --output=results/graph-exp.out
module load gcc/11.2.0

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
            "./${2}" -qtf graphs/${1}/cm${VERTICES}_${DEG}_${i} -s $SEED
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
            "./${2}" -qtf graphs/${1}/er${VERTICES}_${EDGES}_${i} -s $SEED
        done
    done
}

do_topology() {
    for VERTICES in 16 24 32 40 48
    do
        echo "${2} VERTICES: ${VERTICES}"
        for i in {0..5}
        do
            "./${3}" -qtf graphs/${1}/${2}${VERTICES}_*_${i} -s $SEED
        done
    done
}

do_cycle_or_complete() {
    for VERTICES in 16 24 32 40 48
    do
        echo "${2} VERTICES: ${VERTICES}"
        for i in {0..5}
        do
            "./${3}" -qtf graphs/${1}/${2}${VERTICES} -s $SEED
        done
    done
}

# TOPOLOGY
FOLDER="topology_variance"

# SEQ
BIN="digraph16"
RES_PTH="results/topo/seq"
mkdir -p ${RES_PTH}
do_topology $FOLDER "er" $BIN > ${RES_PTH}/er
do_topology $FOLDER "cm" $BIN > ${RES_PTH}/cm
do_cycle_or_complete $FOLDER "c" $BIN > ${RES_PTH}/c
do_cycle_or_complete $FOLDER "k" $BIN > ${RES_PTH}/k

# PAR
BIN="digraph16-PAR"
RES_PTH="results/topo/par"
mkdir -p ${RES_PTH}
do_topology $FOLDER "er" $BIN > ${RES_PTH}/er
do_topology $FOLDER "cm" $BIN > ${RES_PTH}/cm
do_cycle_or_complete $FOLDER "c" $BIN > ${RES_PTH}/c
do_cycle_or_complete $FOLDER "k" $BIN > ${RES_PTH}/k

# EDGE SCALABILITY
FOLDER="edge_scalability"
BIN="digraph16"
RES_PTH="results/edge/seq"
mkdir -p ${RES_PTH}
edge_scale_config $FOLDER $BIN > ${RES_PTH}/cm
edge_scale_er $FOLDER $BIN > ${RES_PTH}/er

# FF COMPARE
FOLDER="ff_comparison"

# SEQ
RES_PTH="results/ff_compare/seq"
mkdir -p ${RES_PTH}
for b in 0 16 32
do
    do_topology $FOLDER "cm" digraph${b} > ${RES_PTH}/digraph${b}
done

# PAR
RES_PTH="results/ff_compare/par"
mkdir -p ${RES_PTH}
for b in 0 16 32
do
    do_topology $FOLDER "cm" digraph${b}-PAR > ${RES_PTH}/digraph${b}
done

# VECTOR PERF
FOLDER="vector_performance"

# SEQ
RES_PTH="results/vecotr_perf/seq"
mkdir -p ${RES_PTH}
do_topology $FOLDER "cm" digraph16 > ${RES_PTH}/vec
do_topology $FOLDER "cm" digraph16-NOVEC > ${RES_PTH}/no_vec

# PAR
RES_PTH="results/vecotr_perf/par"
mkdir -p ${RES_PTH}
do_topology $FOLDER "cm" digraph16-PAR > ${RES_PTH}/vec
do_topology $FOLDER "cm" digraph16-NOVEC-PAR > ${RES_PTH}/no_vec

zip -r $(date '+%Y-%m-%d')_results.zip results
