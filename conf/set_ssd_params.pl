#!/usr/bin/perl

use strict;

# This is optimized for sequential I/O on SSDs.
my $max_block_size = 4096;

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

sub get_main_devname {
	my $dev = $_;
	if ($dev =~ /(sd[a-z]+)[0-9]+/) {
		return $1;
	}
	else {
		return $dev;
	}
}

foreach (@dev_names) {
	my $dev_name = $_;
	my $main_name = get_main_devname($_);
	# sda is the root disk.
	if ($main_name eq "sda") {
		next;
	}

	system("echo noop > /sys/block/${main_name}/queue/scheduler");
	system("echo 2 > /sys/block/${main_name}/queue/rq_affinity");
	system("echo 0 > /sys/block/${main_name}/queue/add_random");
	system("echo $max_block_size > /sys/block/${main_name}/queue/max_sectors_kb");
}
