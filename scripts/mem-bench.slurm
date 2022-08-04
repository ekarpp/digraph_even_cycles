#!/bin/bash

#SBATCH --time=00:05:00
#SBATCH --mem=10G
#SBATCH --partition=short-hsw
#SBATCH --cpus-per-task=24
#SBATCH --ntasks=1
#SBATCH --output=mem-bench.out
module load gcc/11.2.0

hostname
uname -a
cat /etc/*release
g++ --version

make clean
make mem-bench

./mem-bench
