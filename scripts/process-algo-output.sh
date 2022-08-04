#!/bin/bash

REPEATS=6

for f in $(find ./results/*/ -type f -print)
do
    echo $f
    grep -A $(( $REPEATS*3 )) "VERTICES" $f | grep --line-buffered computed \
        | awk 'NR%6 != 1 { printf ("%s %s\n", $4, $7) }' > $f.plot
done

