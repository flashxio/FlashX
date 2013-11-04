#!/usr/bin/perl

use strict;

my $num_args = @ARGV;
if ($num_args < 2) {
	print STDERR "init_test root_conf data_file\n";
	exit 1;
}
my $test_roots = $ARGV[0];
my $data_file = $ARGV[1];
open(FILE, $test_roots) or die "can't open $test_roots: $!";

while (<FILE>) {
	if (/[0-9]+:(.+)/) {
		system("mkdir -p $1");
		system("touch $1/$data_file");
	}
}
