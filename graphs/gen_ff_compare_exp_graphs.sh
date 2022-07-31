#!/bin/bash

REPEATS=6
DEG=5
FOLDER=ff_comparison

mkdir -p $FOLDER

for v in 16 24 32 40 48 56 64
do
    ./config_model.py $v $DEG $REPEATS $FOLDER
done
