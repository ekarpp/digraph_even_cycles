#!/bin/bash

NUM_GRAPHS=6
DEGREE=6

for v in 16 24 32 40 48 56 64
do
    ./c_gen.py $v
    ./k_gen.py $v
    ./erdos_renyi.py $v $(( $v*($v-1) / 2 )) $NUM_GRAPHS
    ./config_model.py $v $DEGREE $NUM_GRAPHS
done

mkdir -p topology_variance
mv config/* topology_variance
mv complete/* topology_variance
mv cycle/* topology_variance
mv erdos_renyi/* topology_variance
