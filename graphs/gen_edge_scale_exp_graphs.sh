#!/bin/bash

REPEATS=6
VERTICES=40
EDGES=$(( $VERTICES * ($VERTICES - 1) ))

while [ $EDGES -gt 1 ]
do
    ./erdos_renyi.py $VERTICES $EDGES $REPEATS
    EDGES=$(echo "$EDGES" | awk '{ print( int($1 / sqrt(2)) ) }')
done

mkdir -p edge_scalability
mv erdos_renyi/* edge_scalability
