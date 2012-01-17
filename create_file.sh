#!/bin/sh

if [ $# -eq 0 ]; then
	echo "specify the filesystem"
	exit 1
fi

umount /mnt/ram
if [ $1 = "ext4" ]; then
	mkfs.ext4 -O extent /dev/ram0
	mount -o dioread_nolock /dev/ram0 /mnt/ram
elif [ $1 = "xfs" ]; then
	mkfs.xfs -f /dev/ram0
	mount /dev/ram0 /mnt/ram
fi

./create_file /mnt/ram/test 4G
chmod 444 /mnt/ram/test
echo 3 > /proc/sys/vm/drop_caches
