#!/bin/bash

#SBATCH --time=1:00:00
#SBATCH --mem=52G
#SBATCH --partition=short-hsw
#SBATCH --cpus-per-task=24
#SBATCH --ntasks=1
#SBATCH --output=results/ext-perf.out
module load gcc/11.2.0

hostname
uname -a
cat /etc/*release
g++ --version

make clean
make objects -j24
make digraph-tests
make extension-perf

make test

ARGS="-s123 -t 30"

for b in 0 16 32
do
    echo "$b bits"
    echo "SEQ"
    ./extension-perf${b} $ARGS
    echo "PAR"
    ./extension-perf${b}-PAR $ARGS
done
