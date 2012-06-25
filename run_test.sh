#!/bin/sh

RAMDISK=/mnt/ram/test
SSD=/dev/sdc	# I assume sdc is an SSD

# test parted global cache
./rand-read /mnt/ram/test option=parted_global pages=1048576 threads=1 workload=RAND_PERMUTE entry_size=4096 cache_size=512M cache_type=associative num_nodes=1
./rand-read /mnt/ram/test option=parted_global pages=1048576 threads=8 workload=RAND_PERMUTE entry_size=4096 cache_size=512M cache_type=associative num_nodes=1

# test global cache
./rand-read $RAMDISK option=global_cache pages=262144 threads=1 workload=RAND_PERMUTE cache_size=512M cache_type=associative entry_size=4096
./rand-read $RAMDISK option=global_cache pages=262144 threads=8 workload=RAND_PERMUTE cache_size=512M cache_type=associative entry_size=4096
./rand-read $SSD option=global_cache pages=40960 threads=1 workload=RAND_PERMUTE cache_size=16M cache_type=associative
