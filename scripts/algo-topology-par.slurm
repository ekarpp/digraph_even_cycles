#!/bin/bash

#SBATCH --time=1:00:00
#SBATCH --mem=1G
#SBATCH --partition=batch-skl
#SBATCH --cpus-per-task=40
#SBATCH --ntasks=1
#SBATCH --output=algo-sqale-par.out

hostname
uname -a
cat /etc/*release
g++ --version

make clean
make digraph-tests16
make digraph16

make test16

REPEATS=6
DEG=6

SEED=123

for v in 16 24 32 40 48 56 64
do
    rm graphs/erdos_renyi/er*
    rm graphs/config/cm*

    echo $v

    # MAKE THESE BEFORE
    cd graphs
    ./k_gen.py $v
    ./c_gen.py $v
    ./config_model.py $v $DEG $REPEATS
    ./erdos_renyi.py $v $(( $v*$v )) $REPEATS
    cd ..

    echo "COMPLETE"
    for (( i=1; i<=$REPEATS; i++ ))
    do
	./digraph16-PAR -qtf graphs/complete/k$v -s $SEED
    done

    echo "CYCLE"
    for (( i=1; i<=$REPEATS; i++ ))
    do
	./digraph16-PAR -qtf graphs/cycle/c$v -s $SEED
    done

    echo "CONFIG"
    for file in graphs/config/*
    do
	./digraph16-PAR -uqtf $file -s $SEED
    done

    echo "ERDOS"
    for file in graphs/erdos_renyi/*
    do
	./digraph16-PAR -uqtf $file -s $SEED
    done
done
