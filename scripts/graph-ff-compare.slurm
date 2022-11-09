#!/bin/bash

# time not accurate
#SBATCH --time=4:00:00
#SBATCH --mem=1G
#SBATCH --partition=short-hsw
#SBATCH --cpus-per-task=24
#SBATCH --ntasks=1
#SBATCH --output=results/graph-exp-ff-compare.out
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

source scripts/helpers

FOLDER="ff_comparison"
RES_PTH="results/ff_compare"
mkdir -p ${RES_PTH}

for b in 0,16,32
do
    do_topology $FOLDER "cm" digraph${b} $SEED > ${RES_PTH}/seq/digraph${b}
done

for b in 0,16,32
do
    do_topology $FOLDER "cm" digraph${b}-PAR $SEED > ${RES_PTH}/par/digraph${b}
done
