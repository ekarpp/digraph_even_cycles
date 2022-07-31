#!/bin/bash

#SBATCH --time=??:00:00
#SBATCH --mem=1G
#SBATCH --partition=batch-skl
#SBATCH --cpus-per-task=40
#SBATCH --ntasks=1
#SBATCH --output=algo-topology-par.out

hostname
uname -a
cat /etc/*release
g++ --version

make clean
make digraph-tests16
make digraph16

make test16

SEED=123

for G in graphs/topology_variance/*
do
    echo ${G##*/}
    ./digraph16-PAR -qtf $G -s $SEED
done
