#!/bin/bash

REPEATS=6
VERTICES=40
DEG=39
EDGES=$(( $VERTICES*($VERTICES-1) ))
FOLDER=edge_scalability

mkdir -p $FOLDER


./k_gen.py ${VERTICES} .

for i in {0..5}
do
    cp k${VERTICES} $FOLDER/cm${VERTICES}_${DEG}_$i
    cp k${VERTICES} $FOLDER/er${VERTICES}_${EDGES}_$i
done

rm k${VERTICES}

DEG=$(echo "$DEG" | awk '{ print( int($1 / (2^(1/3)) ) ) }')
EDGES=$(echo "$EDGES" | awk '{ print( int( $1 / (2^(1/3)) ) ) }')
while [ $DEG -gt 0 ]
do
    echo "$DEG $EDGES"
    ./config_model.py $VERTICES $DEG $REPEATS $FOLDER loops_n_dupes
    ./erdos_renyi.py $VERTICES $EDGES $REPEATS $FOLDER
    DEG=$(echo "$DEG" | awk '{ print( int($1 / (2^(1/3)) ) ) }')
    EDGES=$(echo "$EDGES" | awk '{ print( int( $1 / (2^(1/3)) ) ) }')
done
