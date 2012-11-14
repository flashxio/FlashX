#!/bin/sh

RAMDISK=/mnt/ram/test
SSD=/dev/sdc	# I assume sdc is an SSD

# test parted global cache
#./rand-read /mnt/ram/test option=parted_global pages=1048576 threads=1 workload=RAND_PERMUTE entry_size=4096 cache_size=512M cache_type=associative num_nodes=1
#./rand-read /mnt/ram/test option=parted_global pages=1048576 threads=8 workload=RAND_PERMUTE entry_size=4096 cache_size=512M cache_type=associative num_nodes=1

# test global cache
#./rand-read $RAMDISK option=global_cache pages=262144 threads=1 workload=RAND_PERMUTE cache_size=512M cache_type=associative entry_size=4096
#./rand-read $RAMDISK option=global_cache pages=262144 threads=8 workload=RAND_PERMUTE cache_size=512M cache_type=associative entry_size=4096
#./rand-read $SSD option=global_cache pages=40960 threads=1 workload=RAND_PERMUTE cache_size=16M cache_type=associative

size=4096
npages=10485760
option=aio
pattern=RAND_PERMUTE

node1=0
#numactl -m $node1 -N $node1 ./rand-read /dev/sdb option=$option pages=$npages threads=1 workload=$pattern entry_size=$size num_nodes=1 &
#numactl -m $node1 -N $node1 ./rand-read /dev/sdc option=$option pages=$npages threads=1 workload=$pattern entry_size=$size num_nodes=1 &
#numactl -m $node1 -N $node1 ./rand-read /dev/sdd option=$option pages=$npages threads=1 workload=$pattern entry_size=$size num_nodes=1 &
#numactl -m $node1 -N $node1 ./rand-read /dev/sde option=$option pages=$npages threads=1 workload=$pattern entry_size=$size num_nodes=1 &
#numactl -m $node1 -N $node1 ./rand-read /dev/sdf option=$option pages=$npages threads=1 workload=$pattern entry_size=$size num_nodes=1 &
#numactl -m $node1 -N $node1 ./rand-read /dev/sdp option=$option pages=$npages threads=1 workload=$pattern entry_size=$size num_nodes=1 &

node2=1
#numactl -m $node2 -N $node2 ./rand-read /dev/sdg option=$option pages=$npages threads=1 workload=$pattern entry_size=$size num_nodes=1 &
#numactl -m $node2 -N $node2 ./rand-read /dev/sdh option=$option pages=$npages threads=1 workload=$pattern entry_size=$size num_nodes=1 &
#numactl -m $node2 -N $node2 ./rand-read /dev/sdi option=$option pages=$npages threads=1 workload=$pattern entry_size=$size num_nodes=1 &
#numactl -m $node2 -N $node2 ./rand-read /dev/sdj option=$option pages=$npages threads=1 workload=$pattern entry_size=$size num_nodes=1 &
#numactl -m $node2 -N $node2 ./rand-read /dev/sdk option=$option pages=$npages threads=1 workload=$pattern entry_size=$size num_nodes=1 &
#numactl -m $node2 -N $node2 ./rand-read /dev/sdq option=$option pages=$npages threads=1 workload=$pattern entry_size=$size num_nodes=1 &

node3=1
numactl -m $node3 -N $node3 ./rand-read /dev/sdl option=$option pages=$npages threads=1 workload=$pattern entry_size=$size num_nodes=1 &
numactl -m $node3 -N $node3 ./rand-read /dev/sdm option=$option pages=$npages threads=1 workload=$pattern entry_size=$size num_nodes=1 &
numactl -m $node3 -N $node3 ./rand-read /dev/sdn option=$option pages=$npages threads=1 workload=$pattern entry_size=$size num_nodes=1 &
