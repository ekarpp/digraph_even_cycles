#!/bin/bash

REPEATS=6
VERTICES=40
EDGES=$(( $VERTICES * ($VERTICES - 1) ))
FOLDER=edge_scalability

mkdir -p $FOLDER

while [ $EDGES -gt 1 ]
do
    ./erdos_renyi.py $VERTICES $EDGES $REPEATS $FOLDER
    EDGES=$(echo "$EDGES" | awk '{ print( int($1 / sqrt(2)) ) }')
done
