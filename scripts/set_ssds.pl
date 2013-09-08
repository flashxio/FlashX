#!/usr/bin/perl

use strict;

# this file mounts SSDs to their corresponding directories and sets
# up data files correctly.

# We assign a PCI bus to a NUMA node.
# I hope there are few PCI buses than NUMA nodes.
my $node_id = -1;
my $curr_host = -1;

my @dev_details = `ls -l /sys/block/ | grep sd | grep -v sda`;
my $dev_idx = 1;

open(my $fh, ">", "conf/data_files.txt");
while (@dev_details) {
	my $item = shift(@dev_details);
	if ($item =~ /host([0-9]+).*sd([b-z])$/) {
		my $host = $1;
		my $dev = $2;

		if ($curr_host != $host) {
			$curr_host = $host;
			$node_id++;
		}
		my $data_file = "$node_id:/mnt/ssd${dev_idx}/test_RAID0_16";
		system("echo $data_file > conf/data_file${dev_idx}.txt");
		print $fh "$data_file\n";

		my $dev_file = "/dev/sd${dev}1";
		my $mount_dir = "/mnt/ssd${dev_idx}";
		print "dev: $dev_file, dir: $mount_dir\n" ;
		system("mkdir -p $mount_dir");
		system("mount $dev_file $mount_dir");
		system("chown -R zhengda.zhengda $mount_dir");
		system("echo noop > /sys/block/sd${dev}/queue/scheduler");
		system("cat /sys/block/sd${dev}/queue/scheduler");
		$dev_idx++;
	}
}
close $fh;

$dev_idx--;
print "There are $dev_idx devices\n";

