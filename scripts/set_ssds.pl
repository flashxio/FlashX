#!/bin/sh

use strict;

my @dev_names=("b", "c", "d", "e", "f", "g", "h", "i", "j", "k", "l", "m", "n", "o", "p");
my $file_size = "20G";

for (my $i = 1; $i <= 15; $i++) {
	my $dev_idx = $i - 1;
	my $dev_file = "/dev/sd$dev_names[$dev_idx]1";
	my $dir = "/mnt/ssd${i}";
	print "dev: $dev_file, dir: $dir\n" ;
	system("mkdir -p $dir");
	system("mount $dev_file $dir");
	system("chown zhengda.zhengda $dir/test");
}

