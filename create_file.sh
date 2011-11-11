mkfs.ext4 -O extent /dev/ram0
mount -o dioread_nolock /dev/ram0 /mnt/ram
./create_file /mnt/ram/test 4G
