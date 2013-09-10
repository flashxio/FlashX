#!/usr/bin/perl

# This perl script prints out the backtrace of all threads of a process
# with gdb.

use strict;

my $num_args = @ARGV;
if ($num_args < 2) {
	print STDERR "check_threads.pl program pid\n";
	exit;
}

my $prog = $ARGV[0];
my $pid = $ARGV[1];

system("echo \"info threads\" > /tmp/gdb_cmds.$$.txt");
my @lines = `gdb $prog $pid -x /tmp/gdb_cmds.$$.txt -batch`;
for my $line (@lines) {
	if ($line =~ /^\*?\s+([0-9]+)/) {
		system("echo \"thread $1\" >> /tmp/gdb_cmds_each_thread.$$.txt");
		system("echo backtrace >> /tmp/gdb_cmds_each_thread.$$.txt");
	}
}
system("gdb $prog $pid -x /tmp/gdb_cmds_each_thread.$$.txt -batch");
