#!/bin/sh

dev_files=`ls /dev/nvme*n1`
i=1
rm -f /FlashX/data_files.txt
for file in $dev_files
do
	echo $file
	mkfs.xfs -f $file
	mkdir -p /mnt/ssd$i
	mount $file /mnt/ssd$i
	echo "0:/mnt/ssd$i" >> /FlashX/data_files.txt
	i=$((i+1))
done

jupyter notebook --ip=0.0.0.0 --no-browser --allow-root
