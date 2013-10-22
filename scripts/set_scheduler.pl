#!/usr/bin/perl

use strict;

my @dev_details = `ls -l /sys/block/ | grep sd | grep -v sda`;

my @dev_names;
while (@dev_details) {
	my $item = shift(@dev_details);
	if ($item =~ /sd([b-z])$/) {
		my $dev = $1;
		# get all devices
		push(@dev_names, $dev);
		print "dev: $dev\n";
	}
}

my $scheduler = "noop";

my $num_devs = @dev_names;
for (my $i = 1; $i <= $num_devs; $i++) {
	my $dev_idx = $i - 1;
	my $dev_file = "/dev/sd$dev_names[$dev_idx]";
	system("echo $scheduler > /sys/block/sd$dev_names[$dev_idx]/queue/scheduler");
	system("cat /sys/block/sd$dev_names[$dev_idx]/queue/scheduler");
	system("echo 2 > /sys/block/sd$dev_names[$dev_idx]/queue/rq_affinity");
	system("cat /sys/block/sd$dev_names[$dev_idx]/queue/rq_affinity");
}

