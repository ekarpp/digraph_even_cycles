#!/bin/bash

REPEATS=6
VERTICES=40
DEG=39
FOLDER=edge_scalability

mkdir -p $FOLDER

while [ $DEG -gt 0 ]
do
    ./config_model.py $VERTICES $DEG $REPEATS $FOLDER yes
    DEG=$(echo "$DEG" | awk '{ print( int($1 / sqrt(2)) ) }')
done
