#!/bin/sh

echo "unit-tests"
unit-test/test_open_close

echo "this is the basic test for remote IO."
./test/test_rand_io test/conf/run_remote.txt test

echo "this is the basic test for remote IO with verification enabled."
./test/test_rand_io test/conf/run_remote_virt.txt test1 test2

echo "this is the basic test for remote IO for writes to virtual SSDs."
./test/test_rand_io test/conf/run_remote_virt.txt test1 test2 read_percent=0

echo "this is the basic test for remote IO. Large read."
echo "check the average request size. It should be 128 sectors."
./test/test_rand_io test/conf/run_remote.txt test entry_size=$((4096 * 32)) RAID_block_size=64K

echo "this is the basic test for remote IO read from virtual SSDs with verification enabled. Large read."
./test/test_rand_io test/conf/run_remote_virt.txt test1 test2 entry_size=$((4096 * 32)) RAID_block_size=64K

echo "this is the basic test for remote IO write to virtual SSDs with verification enabled. Large write."
./test/test_rand_io test/conf/run_remote_virt.txt test1 test2 entry_size=$((4096 * 32)) read_percent=0 RAID_block_size=64K

echo "this is the basic test for remote IO. Large read."
echo "check the average request size. It should be 256 sectors. This is to test the code path of splitting a request."
./test/test_rand_io test/conf/run_remote.txt test entry_size=$((4096 * 32))


echo "the basic test for global cached read."
./test/test_rand_io test/conf/run_cache.txt test

echo "the basic test for global cached read on virtual SSDs"
./test/test_rand_io test/conf/run_cache_virt.txt test1 test2

echo "the basic test for global cached write on virtual SSDs"
./test/test_rand_io test/conf/run_cache_virt.txt test1 test2 read_percent=0

echo "the basic test for global cached IO with small reads on virtual SSDs"
./test/test_rand_io test/conf/run_cache_virt.txt test1 test2 pages=$((1024 * 1024)) entry_size=$((128))

echo "the basic test for global cached IO with large reads"
./test/test_rand_io test/conf/run_cache.txt test entry_size=$((4096 * 32)) RAID_block_size=64K

echo "the basic test for global cached IO with large reads on virtual SSDs"
./test/test_rand_io test/conf/run_cache_virt.txt test1 test2 entry_size=$((4096 * 32)) RAID_block_size=64K

echo "the basic test for global cached IO with small writes on virtual SSDs. Large writes and cache flusher are enabled by default."
./test/test_rand_io test/conf/run_cache_virt.txt test1 test2 pages=$((1024 * 1024)) entry_size=$((128)) read_percent=0

echo "the basic test for global cached IO with large writes on virtual SSDs."
./test/test_rand_io test/conf/run_cache_virt.txt test1 test2 read_percent=0 entry_size=$((4096 * 32)) RAID_block_size=64K

echo "the basic test for global cached IO under the TPCC workload."
./test/test_rand_io test/conf/run_cache_tpcc.txt test1 test2

echo "the basic test for global cached IO under the Neo4j workload."
./test/test_rand_io test/conf/run_cache_tpcc.txt test1 test2 workload=test/workload/DijkstraSearch-reconstructed.data


echo "the basic test for parted global cached read on virtual SSDs."
./test/test_rand_io test/conf/run_parted_cache_virt.txt test1 test2

echo "the basic test for parted global cached IO under the TPCC workload on virtual SSDs."
./test/test_rand_io test/conf/run_parted_cache_tpcc.txt test1 test2

echo "the basic test for parted global cached IO under the TPCC workload with few parallel IOs on virtual SSDs."
./test/test_rand_io test/conf/run_parted_cache_tpcc.txt test io_depth=32
