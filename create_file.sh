mkfs.ext4 -O extent /dev/ram0
mount -o dioread_nolock /dev/ram0 /mnt/ram
./create_file /mnt/ram/test 4G
#dd if=/dev/zero of=/mnt/ram/test bs=4096 count=1048576
chmod 444 /mnt/ram/test
