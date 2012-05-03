#!/bin/sh
for ks in 0 1 2 3 4 5 6 7 8 9 10 11 12 13 14 15 16 17 18 19 20 21 22 23 24 25 #26 27 28 29 30
do
    for th in 8
    do
        output/perf_test $th $ks
    done
done


