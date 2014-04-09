#!/bin/bash

app=$1
outfn=$2

sleep 1
top -d 1 -b -p `pgrep $app` > $outfn
