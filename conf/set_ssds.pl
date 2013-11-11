#!/usr/bin/perl

use strict;

# this file mounts SSDs to their corresponding directories and sets
# up data files correctly.

my $num_args = @ARGV;
print "$num_args\n";
if ($num_args < 1) {
	print STDERR "set_ssds dev_list\n";
	exit 1;
}

open(my $dev_fh, "<", $ARGV[0]) or die "open $ARGV[0]: $!";
my @dev_names;
while (<$dev_fh>) {
	chomp($_);
	if ($_ eq "") {
		next;
	}
	push(@dev_names, $_);
	print "ssd: $_\n";
}
my $num_devs = @dev_names;
print "There are $num_devs SSDs\n";

my %devices;
my %host_ids;

foreach (@dev_names) {
	my $item = `ls -l /sys/block/ | grep $_`;
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
	system("echo 2 > /sys/block/sd${dev}/queue/rq_affinity");
	$dev_idx++;
}

close $fh;

$dev_idx--;
print "There are $dev_idx devices\n";

