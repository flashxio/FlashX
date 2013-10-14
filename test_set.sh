#!/bin/sh

echo "this is the basic test for remote IO."
./rand-read test/run_remote.txt test

echo "this is the basic test for remote IO with verification enabled."
./rand-read test/run_remote_virt.txt test

echo "this is the basic test for remote IO for writes to virtual SSDs."
./rand-read test/run_remote_virt.txt test read_percent=0

echo "this is the basic test for remote IO. Large read."
echo "check the average request size. It should be 128 sectors."
./rand-read test/run_remote.txt test entry_size=$((4096 * 32)) RAID_block_size=64K

echo "this is the basic test for remote IO read from virtual SSDs with verification enabled. Large read."
./rand-read test/run_remote_virt.txt test entry_size=$((4096 * 32)) RAID_block_size=64K

echo "this is the basic test for remote IO write to virtual SSDs with verification enabled. Large write."
./rand-read test/run_remote_virt.txt test entry_size=$((4096 * 32)) read_percent=0 RAID_block_size=64K

echo "this is the basic test for remote IO. Large read."
echo "check the average request size. It should be 256 sectors. This is to test the code path of splitting a request."
./rand-read test/run_remote.txt test entry_size=$((4096 * 32))


echo "the basic test for global cached read."
./rand-read test/run_cache.txt test

echo "the basic test for global cached read on virtual SSDs"
./rand-read test/run_cache_virt.txt test

echo "the basic test for global cached write on virtual SSDs"
./rand-read test/run_cache_virt.txt test read_percent=0

echo "the basic test for global cached IO with small reads on virtual SSDs"
./rand-read test/run_cache_virt.txt test pages=$((1024 * 1024)) entry_size=$((128))

echo "the basic test for global cached IO with large reads"
./rand-read test/run_cache.txt test entry_size=$((4096 * 32)) RAID_block_size=64K

echo "the basic test for global cached IO with large reads on virtual SSDs"
./rand-read test/run_cache_virt.txt test entry_size=$((4096 * 32)) RAID_block_size=64K

echo "the basic test for global cached IO with small writes on virtual SSDs. Large writes and cache flusher are enabled by default."
./rand-read test/run_cache_virt.txt test pages=$((1024 * 1024)) entry_size=$((128)) read_percent=0

echo "the basic test for global cached IO with large writes on virtual SSDs."
./rand-read test/run_cache_virt.txt test read_percent=0 entry_size=$((4096 * 32)) RAID_block_size=64K
