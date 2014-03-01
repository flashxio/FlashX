#!/bin/sh

# unit tests

echo "unit-tests"
unit-test/test_open_close


# correctness test with virtual IO.

echo "this is the basic test for remote IO with verification enabled."
./test/test_rand_io test/conf/run_remote_virt.txt test1 test2

echo "this is the basic test for remote IO for writes to virtual SSDs."
./test/test_rand_io test/conf/run_remote_virt.txt test1 test2 read_percent=0

echo "this is the basic test for remote IO read from virtual SSDs with verification enabled. Large read."
./test/test_rand_io test/conf/run_remote_virt.txt test1 test2 entry_size=$((4096 * 32)) RAID_block_size=64K

echo "this is the basic test for remote IO write to virtual SSDs with verification enabled. Large write."
./test/test_rand_io test/conf/run_remote_virt.txt test1 test2 entry_size=$((4096 * 32)) read_percent=0 RAID_block_size=64K


# max speed: read 320.80MB/s, write 0.00MB/s
echo "the basic test for global cached read on virtual SSDs"
./test/test_rand_io test/conf/run_cache_virt.txt test1 test2
./test/test_rand_io test/conf/run_cache_virt.txt test1 test2 user_compute=

# max speed: read 0.00MB/s, write 254.75MB/s
echo "the basic test for global cached write on virtual SSDs"
./test/test_rand_io test/conf/run_cache_virt.txt test1 test2 read_percent=0
./test/test_rand_io test/conf/run_cache_virt.txt test1 test2 read_percent=0 user_compute=

# max speed: read 307.26MB/s, write 0.00MB/s
echo "the basic test for global cached IO with large reads on virtual SSDs"
./test/test_rand_io test/conf/run_cache_virt.txt test1 test2 entry_size=$((4096 * 32)) RAID_block_size=64K
./test/test_rand_io test/conf/run_cache_virt.txt test1 test2 entry_size=$((4096 * 32)) RAID_block_size=64K user_compute=

# max speed: read 0.00MB/s, write 556.09MB/s
echo "the basic test for global cached IO with large writes on virtual SSDs."
./test/test_rand_io test/conf/run_cache_virt.txt test1 test2 read_percent=0 entry_size=$((4096 * 32)) RAID_block_size=64K
./test/test_rand_io test/conf/run_cache_virt.txt test1 test2 read_percent=0 entry_size=$((4096 * 32)) RAID_block_size=64K user_compute=

# max speed: read 192.80MB/s, write 144.09MB/s
echo "the basic test for global cached IO under the TPCC workload."
./test/test_rand_io test/conf/run_cache_tpcc.txt test1 test2
./test/test_rand_io test/conf/run_cache_tpcc.txt test1 test2 user_compute=

# max speed: read 259.47MB/s, write 0.00MB/s
echo "the basic test for global cached IO under the Neo4j workload."
./test/test_rand_io test/conf/run_cache_tpcc.txt test1 test2 workload=test/workload/DijkstraSearch-reconstructed.data
./test/test_rand_io test/conf/run_cache_tpcc.txt test1 test2 workload=test/workload/DijkstraSearch-reconstructed.data user_compute=

# max speed: read 285.50MB/s, write 0.00MB/s
# cache hits: 254,499,141
# total time: 793.357666 seconds
echo "the basic test for global cached IO under my own triangle counting workload"
./test/test_rand_io test/conf/run_cache_tpcc.txt test1 test2 workload=test/workload/triangle-counting.data
./test/test_rand_io test/conf/run_cache_tpcc.txt test1 test2 workload=test/workload/triangle-counting.data user_compute=


# max speed: read 235.30MB/s, write 0.00MB/s
echo "the basic test for parted global cached read on virtual SSDs."
./test/test_rand_io test/conf/run_parted_cache_virt.txt test1 test2
./test/test_rand_io test/conf/run_parted_cache_virt.txt test1 test2 user_compute=

# max speed: read 171.32MB/s, write 85.96MB/s
echo "the basic test for parted global cached IO under the TPCC workload on virtual SSDs."
./test/test_rand_io test/conf/run_parted_cache_tpcc.txt test1 test2

echo "the basic test for parted global cached IO under the TPCC workload with few parallel IOs on virtual SSDs."
./test/test_rand_io test/conf/run_parted_cache_tpcc.txt test io_depth=32


# performance test on the real SSDs.

# total time: 14.2s to read 40GB.
echo "the basic test for remote read."
./test/test_rand_io test/conf/run_remote.txt test

# total time: 8.1s to read 40GB.
echo "this is the basic test for remote IO. Large read."
echo "check the average request size. It should be 128 sectors."
./test/test_rand_io test/conf/run_remote.txt test entry_size=$((4096 * 32)) RAID_block_size=64K

# total time: 13s to read 40GB
echo "the basic test for global cached read."
./test/test_rand_io test/conf/run_cache.txt test

# total time: 7.8s to read 40GB.
echo "the basic test for global cached IO with large reads"
./test/test_rand_io test/conf/run_cache.txt test entry_size=$((4096 * 32)) RAID_block_size=64K

echo "the basic test for global cached IO under the TPCC workload."
./test/test_rand_io test/conf/run_cache_real.txt test workload=test/workload/mysqld_tpcc_datafile.data

# total time: 10.8s
# cache hits: 19296954/26326945
echo "the basic test for global cached IO under the Neo4j workload."
./test/test_rand_io test/conf/run_cache_real.txt test workload=test/workload/DijkstraSearch-reconstructed.data cache_size=512M

# total time: 273s
# cache hits: 643938618/829040171
echo "the basic test for global cached IO under my own triangle counting workload"
./test/test_rand_io test/conf/run_cache_real.txt test workload=test/workload/triangle-counting.data

# total time: 268s
./test/test_rand_io test/conf/run_cache_real.txt test workload=test/workload/triangle-counting.data user_compute=

./test/test_rand_io test/conf/run_parted_cache_real.txt test workload=test/workload/triangle-counting.data user_compute=

apps/graph-bfs/graph-bfs test/conf/run_graph.txt web-graph-v3 web-graph-index-v3 0
apps/sssp/sssp test/conf/run_graph.txt web-graph-v3 web-graph-index-v3 0
apps/triangle-counting/triangle-counting test/conf/run_graph.txt soc-LiveJournal-v3 soc-LiveJournal-index-v3
apps/sstsg/scan-statistics test/conf/run_graph.txt bitcoin-v3 bitcoin-index-v3 50 5
apps/scan-statistics/scan-statistics test/conf/run_graph.txt web-graph-v3 web-graph-index-v3 -c "max_processing_vertices=3"
