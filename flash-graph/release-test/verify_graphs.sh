#!/bin/sh

# verify generating large graphs.

../tools/el2al -v -w -T 32 -d /mnt/ram4/twitter.adj /mnt/ram4/twitter.index /mnt/nfs1/graph-data/twitter_rv.net.gz
rm /mnt/ram4/*
../tools/el2al -u -v -w -T 32 -d /mnt/ram4/friendster.adj /mnt/ram4/friendster.index /mnt/nfs1/graph-data/com-friendster.ungraph.txt.gz
rm /mnt/ram4/*
