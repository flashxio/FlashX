#!/usr/bin/perl

use strict;
use POSIX;
use Time::HiRes qw(usleep nanosleep);

my $enable_log = 1;
my $num_args = @ARGV;
if ($num_args < 1) {
	print STDERR "monitor_proc cmd\n";
	exit;
}
my $log_file = "monitor.log";
if ($enable_log) {
	open LOG, ">", $log_file or die "can't open $log_file: $!";
}
my $cmd = $ARGV[0];
my $pid = run($cmd);
while (1) {
	print "========================================\n";
	my @tids = collect_info($pid);
	# the process probably has exit.
	if (@tids == 0) {
		last;
	}
	waitpid $pid, WNOHANG;
	usleep(300000);
}
if ($enable_log) {
	close LOG;
}

sub parse_numa_maps {
	my @lines = @{$_[0]};
	my @mem_sizes = (0, 0, 0, 0);
	for (my $i = 0; $i < @lines; $i++) {
		my @fields = split / /, $lines[$i];
		for (my $j = 0; $j < @fields; $j++) {
			if ($fields[$j] =~ /N([0-9]+)=([0-9]+)/) {
				$mem_sizes[$1] += $2;
			}
		}
	}
	for (my $i = 0; $i < @mem_sizes; $i++) {
		print "node$i: $mem_sizes[$i]\n"
	}
}

sub collect_info {
	my $pid = $_[0];
	my $numa_maps = `cat /proc/$pid/numa_maps`;
	if ($enable_log) {
		print LOG $numa_maps;
		print LOG `cat /proc/$pid/maps`;
	}
	my @lines = split /\n/, $numa_maps;
	parse_numa_maps(\@lines);

	my @tids = collect_tids($pid);
	for (my $i = 0; $i < @tids; $i++) {
		my @fields = split(/ /, `cat /proc/$pid/task/$tids[$i]/stat`);
		print "thread $tids[$i] runs on core $fields[38]\n";
	}
	return @tids;
}

sub collect_tids {
	my $pid = $_[0];
	my @tids;
	opendir(D, "/proc/$pid/task") || return ();
	while (my $f = readdir(D)) {
		if (isdigit($f)) {
			push(@tids, $f);
		}
	}
	closedir(D);
	return @tids;
}

sub run {
	my $cmd = $_[0];
	my $pid = fork();
	if ($pid == 0) {		# child
		exec($cmd) or die "can't execute: $!";
		exit;
	}
	elsif ($pid > 0) {
		return $pid;
	}
	else {
		die "can't create a process: $!";
	}
}
