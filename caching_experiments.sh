#!/bin/bash

if [ "$#" -ne 2 ]; then
  echo "usage: caching_experiments.sh graph-file graph-index-file"
  exit -1
fi

conf="test/conf/run_graph.txt"
outdir=CacheTests

graph=$1
graph_index=$2
#graph=page-graph-v3
#graph_index=page-graph-index-v3

for CACHE_SIZE in 4 16 64 256
do
  echo "Running bfs with cache size = $CACHE_SIZE... "
  ./sleep_top.sh graph-bfs $outdir/bfs_"$CACHE_SIZE"gb_cache.out &
  apps/graph-bfs/graph-bfs $conf $graph $graph_index 0 -c "cache_size="$CACHE_SIZE"G"
  kill `pgrep sleep_top` 

  echo "Running triangle counting with cache size = $CACHE_SIZE... "
  ./sleep_top.sh triangle-counting $outdir/tri_"$CACHE_SIZE"gb_cache.out &
  apps/triangle-counting/triangle-counting $conf $graph $graph_index -c "cache_size="$CACHE_SIZE"G"
  kill `pgrep sleep_top` 

  echo "Running scan stat with cache size = $CACHE_SIZE... "
  ./sleep_top.sh topK-scan $outdir/scan_"$CACHE_SIZE"gb_cache.out &
  apps/scan-statistics/topK-scan $conf $graph $graph_index 0 -c "cache_size="$CACHE_SIZE"G max_processing_vertices=2"
  kill `pgrep sleep_top` 

  echo "Running wcc with cache size = $CACHE_SIZE... "
  ./sleep_top.sh wcc $outdir/wcc_"$CACHE_SIZE"gb_cache.out &
  apps/weakly-connected-components/wcc $conf $graph $graph_index 0 -s 1000000 -o $outdir/tmp.txt -c "cache_size="$CACHE_SIZE"G"

  kill `pgrep sleep_top` 
  echo "Running page-rank2 with cache size = $CACHE_SIZE... "
  ./sleep_top.sh page-rank2 $outdir/pgrnk_"$CACHE_SIZE"gb_top_out.txt &
  apps/page-rank/page-rank2 $conf $graph $graph_index 0.85 -c "cache_size="$CACHE_SIZE"G" -i 30
  kill `pgrep sleep_top` 

done
