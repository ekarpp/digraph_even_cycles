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

make test16
make digraph16

for vertices in 16 24 32 40 48
do
    echo "VECTORIZED vertices: ${vertices}"
    for i in {0..5}
    do
        graph=graphs/ff_comparison/cm${vertices}_*_${i}
        ./digraph16 -qtf $graph -s $SEED
    done

    echo "nonvecs vertices: ${vertices}"
    for i in {0..5}
    do
        graph=graphs/ff_comparison/cm${vertices}_*_${i}
        ./digraph16-NOVEC -qtf $graph -s $SEED
    done

done
