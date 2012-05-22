#!/bin/sh

# test parted global cache
./rand-read /mnt/ram/test option=parted_global pages=1048576 threads=1 workload=RAND_PERMUTE entry_size=4096 cache_size=512M cache_type=associative num_nodes=1
./rand-read /mnt/ram/test option=parted_global pages=1048576 threads=8 workload=RAND_PERMUTE entry_size=4096 cache_size=512M cache_type=associative num_nodes=1

# test global cache
./rand-read /mnt/ram/test option=global_cache pages=262144 threads=1 workload=RAND_PERMUTE cache_size=512M cache_type=associative entry_size=4096
./rand-read /mnt/ram/test option=global_cache pages=262144 threads=8 workload=RAND_PERMUTE cache_size=512M cache_type=associative entry_size=4096

