#!/bin/sh

size="40G"

if [ $# > 1 ]; then
	size=$1
fi
echo "file size is $size"

mkdir -p /mnt/ram0
mkdir -p /mnt/ram1
mkdir -p /mnt/ram2
mkdir -p /mnt/ram3

mkfs.xfs /dev/ram0
mkfs.xfs /dev/ram1
mkfs.xfs /dev/ram2
mkfs.xfs /dev/ram3

mount /dev/ram0 /mnt/ram0
mount /dev/ram1 /mnt/ram1
mount /dev/ram2 /mnt/ram2
mount /dev/ram3 /mnt/ram3

numactl -m 0 -N 0 tools/create_file /mnt/ram0/test $size
numactl -m 1 -N 1 tools/create_file /mnt/ram1/test $size
numactl -m 2 -N 2 tools/create_file /mnt/ram2/test $size
numactl -m 3 -N 3 tools/create_file /mnt/ram3/test $size

chown -R zhengda.zhengda /mnt/ram0
chown -R zhengda.zhengda /mnt/ram1
chown -R zhengda.zhengda /mnt/ram2
chown -R zhengda.zhengda /mnt/ram3
