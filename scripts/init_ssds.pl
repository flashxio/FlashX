#!/usr/bin/perl

use strict;

my @dev_names=("b", "c", "d", "e", "f", "g", "h", "i", "j", "k", "l", "m", "n", "o", "p", "q");
my $file_size = "40G";

for (my $i = 1; $i <= 16; $i++) {
	my $dev_idx = $i - 1;
	my $dev_file = "/dev/sd$dev_names[$dev_idx]1";
	my $dir = "/mnt/ssd${i}";
	print "dev: $dev_file, dir: $dir\n" ;
	system("mkfs.xfs $dev_file");
	system("mkdir -p $dir");
	system("mount $dev_file $dir");
}

system("tools/create_files data_files_RAID0.txt 640G RAID0 16");

for (my $i = 1; $i <= 16; $i++) {
	my $dev_idx = $i - 1;
	my $dev_file = "/dev/sd$dev_names[$dev_idx]1";
	my $dir = "/mnt/ssd${i}";
	system("chown -R zhengda.zhengda $dir");
	system("echo noop > /sys/block/sd$dev_names[$dev_idx]/queue/scheduler");
	system("cat /sys/block/sd$dev_names[$dev_idx]/queue/scheduler");
}
