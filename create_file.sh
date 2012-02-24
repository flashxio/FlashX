#!/bin/sh

if [ $# -eq 0 ]; then
	echo "create_file filesystem size ram node_id"
	exit 1
fi

size=$2
ram=$3
node_id=$4

umount /mnt/$ram
mkdir -p /mnt/$ram
if [ $1 = "ext4" ]; then
	mkfs.ext4 -O extent /dev/$ram
	mount -o dioread_nolock /dev/$ram /mnt/$ram
	chmod 444 /mnt/$ram/test
elif [ $1 = "xfs" ]; then
	mkfs.xfs -f /dev/$ram
	mount /dev/$ram /mnt/$ram
elif [ $1 = "tmpfs" ]; then
	mount -t tmpfs -o size=5G tmpfs /mnt/tmp
fi

./create_file /mnt/$ram/test $size $node_id
echo 3 > /proc/sys/vm/drop_caches
