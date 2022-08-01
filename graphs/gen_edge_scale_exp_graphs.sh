#!/bin/bash

REPEATS=6
VERTICES=40
DEG=39
FOLDER=edge_scalability

mkdir -p $FOLDER


./k_gen.py ${VERTICES} .

for i in {0..5}
do
    cp k${VERTICES} $FOLDER/cm${VERTICES}_${DEG}_$i
done

rm k${VERTICES}

DEG=$(echo "$DEG" | awk '{ print( int($1 / sqrt(sqrt(2))) ) }')

while [ $DEG -gt 0 ]
do
    echo $DEG
    ./config_model.py $VERTICES $DEG $REPEATS $FOLDER yes
    DEG=$(echo "$DEG" | awk '{ print( int($1 / sqrt(sqrt(2))) ) }')
done
