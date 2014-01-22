#!/usr/bin/perl

use strict;

my $num_args = @ARGV;
if ($num_args < 1) {
	print STDERR "show_cpu.pl prog_name\n";
	exit;
}
my $prog_name = $ARGV[0];
my @pids = `ps -aux 2>/dev/null | awk '{print \$2, \$11}' | grep $prog_name | awk '{print \$1}'`;
my $num_pids = @pids;
if ($num_pids != 1) {
	print STDERR "There are $num_pids processes named after $prog_name\n";
	foreach(@pids) {
		print $_, "\n";
	}
	exit;
}

system("top -b -p $pids[0] -d 1");
