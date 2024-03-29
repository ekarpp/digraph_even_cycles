#!/bin/bash

NUM_GRAPHS=6
DEGREE=5
FOLDER=topology_variance

mkdir -p $FOLDER

for v in 16 24 32 40 48 56 64
do
    ./c_gen.py $v $FOLDER
    ./k_gen.py $v $FOLDER

    ./erdos_renyi.py $v $(( $v*$DEGREE )) $NUM_GRAPHS $FOLDER
    ./undir_config_model.py $v $DEGREE $NUM_GRAPHS $FOLDER
done
