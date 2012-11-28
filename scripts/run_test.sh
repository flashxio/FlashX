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

size=$((4096 * 1))
npages=$((10485760 * 1))
option=aio
nthreads=1
pattern=RAND_PERMUTE

node1=0
#numactl -m $node1 -N $node1 ./rand-read /dev/sdb option=$option pages=$npages threads=$nthreads workload=$pattern entry_size=$size num_nodes=1 &
#numactl -m $node1 -N $node1 ./rand-read /dev/sdc option=$option pages=$npages threads=$nthreads workload=$pattern entry_size=$size num_nodes=1 &
#numactl -m $node1 -N $node1 ./rand-read /dev/sdd option=$option pages=$npages threads=$nthreads workload=$pattern entry_size=$size num_nodes=1 &
#numactl -m $node1 -N $node1 ./rand-read /dev/sde option=$option pages=$npages threads=$nthreads workload=$pattern entry_size=$size num_nodes=1 &
#numactl -m $node1 -N $node1 ./rand-read /dev/sdf option=$option pages=$npages threads=$nthreads workload=$pattern entry_size=$size num_nodes=1 &

node2=1
numactl -m $node2 -N $node2 ./rand-read /dev/sdc option=$option pages=$npages threads=$nthreads workload=$pattern entry_size=$size num_nodes=1 &
numactl -m $node2 -N $node2 ./rand-read /dev/sdd option=$option pages=$npages threads=$nthreads workload=$pattern entry_size=$size num_nodes=1 &
#numactl -m $node2 -N $node2 ./rand-read /dev/sde option=$option pages=$npages threads=$nthreads workload=$pattern entry_size=$size num_nodes=1 &
#numactl -m $node2 -N $node2 ./rand-read /dev/sdf option=$option pages=$npages threads=$nthreads workload=$pattern entry_size=$size num_nodes=1 &
#numactl -m $node2 -N $node2 ./rand-read /dev/sdg option=$option pages=$npages threads=$nthreads workload=$pattern entry_size=$size num_nodes=1 &
#numactl -m $node2 -N $node2 ./rand-read /dev/sdh option=$option pages=$npages threads=$nthreads workload=$pattern entry_size=$size num_nodes=1 &
#numactl -m $node2 -N $node2 ./rand-read /dev/sdi option=$option pages=$npages threads=$nthreads workload=$pattern entry_size=$size num_nodes=1 &

node3=1
#numactl -m $node3 -N $node3 ./rand-read /dev/sdj option=$option pages=$npages threads=$nthreads workload=$pattern entry_size=$size num_nodes=1 &
#numactl -m $node3 -N $node3 ./rand-read /dev/sdk option=$option pages=$npages threads=$nthreads workload=$pattern entry_size=$size num_nodes=1 &
#numactl -m $node3 -N $node3 ./rand-read /dev/sdl option=$option pages=$npages threads=$nthreads workload=$pattern entry_size=$size num_nodes=1 &
#numactl -m $node3 -N $node3 ./rand-read /dev/sdm option=$option pages=$npages threads=$nthreads workload=$pattern entry_size=$size num_nodes=1 &
#numactl -m $node3 -N $node3 ./rand-read /dev/sdn option=$option pages=$npages threads=$nthreads workload=$pattern entry_size=$size num_nodes=1 &
#numactl -m $node3 -N $node3 ./rand-read /dev/sdo option=$option pages=$npages threads=$nthreads workload=$pattern entry_size=$size num_nodes=1 &
#numactl -m $node3 -N $node3 ./rand-read /dev/sdp option=$option pages=$npages threads=$nthreads workload=$pattern entry_size=$size num_nodes=1 &

node2=1
nthreads=$((14 * 4))
option=aio
#numactl -m $node2 -N $node2 ./rand-read /dev/sdb /dev/sdd /dev/sde /dev/sdg /dev/sdh /dev/sdi /dev/sdj /dev/sdk /dev/sdl /dev/sdm /dev/sdn /dev/sdo /dev/sdp /dev/sdq option=$option pages=$npages threads=$nthreads workload=$pattern entry_size=$size num_nodes=1 &

#numactl -m 1 -N 1 ./rand-read /mnt/ssd2/test option=remote pages=1048576 threads=4 workload=RAND_PERMUTE entry_size=4096 num_nodes=1
#numactl -m 1 -N 1 /home/zhengda/bin/perf record -g ./rand-read /mnt/ssd2/test /mnt/ssd3/test /mnt/ssd4/test /mnt/ssd5/test option=global_cache pages=10485760 threads=4 workload=RAND_PERMUTE entry_size=4096 num_nodes=1 cache_type=associative
#numactl -m 1 -N 1 ./rand-read /mnt/ssd2/test /mnt/ssd3/test /mnt/ssd4/test /mnt/ssd5/test /mnt/ssd6/test /mnt/ssd7/test /mnt/ssd8/test /mnt/ssd9/test /mnt/ssd10/test /mnt/ssd11/test /mnt/ssd12/test /mnt/ssd13/test /mnt/ssd14/test /mnt/ssd15/test option=remote pages=5242880 threads=8 workload=RAND_PERMUTE entry_size=4096 num_nodes=1 verify_content=
