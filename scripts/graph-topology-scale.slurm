#!/bin/bash

# time not accurate
#SBATCH --time=4:00:00
#SBATCH --mem=1G
#SBATCH --partition=short-hsw
#SBATCH --cpus-per-task=24
#SBATCH --ntasks=1
#SBATCH --output=results/graph-exp-topology.out
module load gcc/11.2.0

hostname
uname -a
cat /etc/*release
g++ --version

make clean
make objects -j24
make digraph16
make test16

SEED=123

source scripts/graph_helpers

# TOPOLOGY
FOLDER="topology_variance"

# SEQ
BIN="digraph16"
RES_PTH="results/topo/seq"
mkdir -p ${RES_PTH}
do_topology $FOLDER "er" $BIN $SEED > ${RES_PTH}/er
do_topology $FOLDER "cm" $BIN $SEED > ${RES_PTH}/cm
do_cycle_or_complete $FOLDER "c" $BIN $SEED > ${RES_PTH}/c
do_cycle_or_complete $FOLDER "k" $BIN $SEED > ${RES_PTH}/k

# PAR
BIN="digraph16-PAR"
RES_PTH="results/topo/par"
mkdir -p ${RES_PTH}
do_topology $FOLDER "er" $BIN $SEED > ${RES_PTH}/er
do_topology $FOLDER "cm" $BIN $SEED > ${RES_PTH}/cm
do_cycle_or_complete $FOLDER "c" $BIN $SEED > ${RES_PTH}/c
do_cycle_or_complete $FOLDER "k" $BIN $SEED > ${RES_PTH}/k
