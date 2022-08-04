#!/bin/bash

#SBATCH --time=??:00:00
#SBATCH --mem=1G
#SBATCH --partition=batch-skl
#SBATCH --cpus-per-task=1
#SBATCH --ntasks=1
#SBATCH --output=ff_comparison.out
module load gcc/11.2.0

hostname
uname -a
cat /etc/*release
g++ --version

make objects
make digraph
make test

SEED=123

for bits in 0 16 32
do
    echo "bits: ${bits}"
    for vertices in 16 24 32 40 48
    do
        echo "vertices: ${vertices}"
        for i in {0..5}
        do
            graph=graphs/ff_comparison/cm${vertices}_*_${i}
            ./digraph${b} -qtf $graph -s $SEED
        done

    done
done
