#!/bin/bash

#SBATCH --time=??:00:00
#SBATCH --mem=1G
#SBATCH --partition=batch-skl
#SBATCH --cpus-per-task=1
#SBATCH --ntasks=1
#SBATCH --output=vector-perf.out
module load gcc/11.2.0

hostname
uname -a
cat /etc/*release
g++ --version

make objects
make digraph
make test

for G in graphs/ff_comparison/*
do
        echo $G
done
