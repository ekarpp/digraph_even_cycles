#!/bin/bash

REPEATS=6
VERTICES=40
DEG=39
FOLDER=edge_scalability

mkdir -p $FOLDER

DEG=$(echo "$DEG" | awk '{ print( int($1 / sqrt(2)) ) }')

while [ $DEG -gt 0 ]
do
    echo $DEG
    ./config_model.py $VERTICES $DEG $REPEATS $FOLDER dups
    DEG=$(echo "$DEG" | awk '{ print( int($1 / sqrt(2)) ) }')
done
