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

sub get_main_devname {
	my $dev = $_;
	if ($dev =~ /(xvd[a-z]+)[0-9]+/) {
		return $1;
	}
	else {
		return $dev;
	}
}

my $dev_idx = 1;
open(my $fh, ">", "conf/data_files.txt");
for (@dev_names) {
	my $dev = $_;
	my $main_name = get_main_devname($dev);
	my $data_file = "0:/mnt/ssd${dev_idx}/";
	print $fh "$data_file\n";

	my $dev_file = "/dev/${dev}";
	my $mount_dir = "/mnt/ssd${dev_idx}";
	print "dev: $dev_file, dir: $mount_dir\n" ;
	system("mkdir -p $mount_dir");
	system("mkfs.xfs /dev/$dev");
	system("mount $dev_file $mount_dir");
	system("chown -R ubuntu $mount_dir");
	system("echo noop > /sys/block/${main_name}/queue/scheduler");
	system("cat /sys/block/${main_name}/queue/scheduler");
	system("echo 2 > /sys/block/${main_name}/queue/rq_affinity");
	system("echo 0 > /sys/block/${main_name}/queue/add_random");
	$dev_idx++;
}

close $fh;

$dev_idx--;
print "There are $dev_idx devices\n";
