#!/usr/bin/perl

use strict;

# this file mounts SSDs to their corresponding directories and sets
# up data files correctly.

my @dev_details = `ls -l /sys/block/ | grep sd | grep -v sda`;

my %devices;
my %host_ids;

while (@dev_details) {
	my $item = shift(@dev_details);
	if ($item =~ /host([0-9]+).*sd([b-z])$/) {
		my $host = $1;
		my $dev = $2;

		# get all devices
		$devices{$dev} = $host;
		# get all host Ids.
		$host_ids{$host} = 1;
	}
}

# We assign a PCI bus to a NUMA node.
# I hope there are few PCI buses than NUMA nodes.
my $node_id = 0;
for (sort {$a <=> $b} keys %host_ids) {
	$host_ids{$_} = $node_id;
	print "host $_ is assigned to node $node_id\n";
	$node_id++;
}

my $dev_idx = 1;
open(my $fh, ">", "conf/data_files.txt");
for (sort keys %devices) {
	my $dev = $_;
	my $host = $devices{$dev};
	my $node_id = $host_ids{$host};

	my $data_file = "$node_id:/mnt/ssd${dev_idx}/";
	print $fh "$data_file\n";

	my $dev_file = "/dev/sd${dev}";
	my $mount_dir = "/mnt/ssd${dev_idx}";
	print "dev: $dev_file, dir: $mount_dir, on node $node_id\n" ;
	system("mkdir -p $mount_dir");
	system("mount $dev_file $mount_dir");
	system("chown -R zhengda.zhengda $mount_dir");
	system("echo noop > /sys/block/sd${dev}/queue/scheduler");
	system("cat /sys/block/sd${dev}/queue/scheduler");
	$dev_idx++;
}

close $fh;

$dev_idx--;
print "There are $dev_idx devices\n";

