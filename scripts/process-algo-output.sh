#!/bin/bash

REPEATS=6

for f in $(find ./results/*/ -type f -print)
do
    echo $f
    # everything except edge
    grep -A $(( $REPEATS*3 )) "VERTICES" $f | grep --line-buffered computed \
        | awk 'NR%6 != 1 { printf ("%s %s\n", $4, $7) }' > $f.plot
    # edge
    grep -A $(( $REPEATS*3 )) "degree" $f | awk \
        'NR%19 == 1 { edges=$3*40 }
         NR%19 == 0 || (NR%19 > 6 && NR%19%3 == 1) { printf ("%s %s\n", edges, $7) }' >> $f.plot
done
