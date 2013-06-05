#!/usr/bin/perl

use strict;

my $num_args = @ARGV;

my $iterations = 100000000;
my $interval = 1;	# 1 second
if ($num_args < 0) {
	print "show_cpu.pl node_id\n";
	exit;
}

my $node_id = $ARGV[0];

my @cores;
for (my $i = 0; $i < 8; $i++) {
	push(@cores, $node_id + 4 * $i);
}

my $str = join(',', @cores);
print $str, "\n";

system("mpstat -P $str $interval");
