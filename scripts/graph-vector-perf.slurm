#!/bin/bash

# time not accurate
#SBATCH --time=4:00:00
#SBATCH --mem=1G
#SBATCH --partition=short-hsw
#SBATCH --cpus-per-task=24
#SBATCH --ntasks=1
#SBATCH --output=results/graph-exp-vector.out
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

# VECTOR PERF
FOLDER="vector_performance"
RES_PTH="results/vecotr_perf/seq"
mkdir -p ${RES_PTH}
do_topology $FOLDER "cm" digraph16 $SEED > ${RES_PTH}/vec
do_topology $FOLDER "cm" digraph16-NOVEC $SEED > ${RES_PTH}/no_vec

RES_PTH="results/vecotr_perf/par"
mkdir -p ${RES_PTH}
do_topology $FOLDER "cm" digraph16-PAR $SEED > ${RES_PTH}/vec
do_topology $FOLDER "cm" digraph16-NOVEC-PAR $SEED > ${RES_PTH}/no_vec
