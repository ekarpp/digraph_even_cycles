#!/bin/bash

# time not accurate
#SBATCH --time=4:00:00
#SBATCH --mem=1G
#SBATCH --partition=short-hsw
#SBATCH --cpus-per-task=24
#SBATCH --ntasks=1
#SBATCH --output=results/graph-exp-edge-scale.out
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

# EDGE SCALABILITY
BIN="digraph16"
FOLDER="edge_scalability"
RES_PTH="results/edge/seq"
mkdir -p ${RES_PTH}
edge_scale_config $FOLDER $BIN $SEED > ${RES_PTH}/cm
edge_scale_er $FOLDER $BIN $SEED > ${RES_PTH}/er
