#!/usr/bin/perl

use strict;

my @dev_names=("b", "c", "d", "e", "f", "g", "h", "i", "j", "k", "l", "m", "n", "o", "p", "q");
my $scheduler = "noop";

for (my $i = 1; $i <= 16; $i++) {
	my $dev_idx = $i - 1;
	my $dev_file = "/dev/sd$dev_names[$dev_idx]1";
	system("echo $scheduler > /sys/block/sd$dev_names[$dev_idx]/queue/scheduler");
	system("cat /sys/block/sd$dev_names[$dev_idx]/queue/scheduler");
	system("echo 2 > /sys/block/sd$dev_names[$dev_idx]/queue/rq_affinity");
	system("cat /sys/block/sd$dev_names[$dev_idx]/queue/rq_affinity");
}

