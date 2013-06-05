#!/usr/bin/perl

use strict;
use POSIX;

my $num_args = @ARGV;
if ($num_args < 1) {
	print STDERR "proc_stat pid\n";
	exit;
}
my $pid = $ARGV[0];
my @fields = split(/ /, `cat /proc/$pid/stat`);
print_fields(\@fields);

print "The threads in process $pid:\n";
my @tids = collect_tids($pid);
for (my $i = 0; $i < @tids; $i++) {
	@fields = split(/ /, `cat /proc/$pid/task/$tids[$i]/stat`);
	print "=====================thread $tids[$i]===================\n";
	print_fields(\@fields);
}

sub collect_tids {
	my $pid = $_[0];
	my @tids;
	opendir(D, "/proc/$pid/task") || die "Can't opedir: $!\n";
	while (my $f = readdir(D)) {
		if (isdigit($f)) {
			push(@tids, $f);
		}
	}
	closedir(D);
	return @tids;
}

sub print_fields {
	my @fields = @{$_[0]};
	print "1.pid: $fields[0]\n";
	print "2.executable filename: $fields[1]\n";
	print "3.process state: $fields[2]\n";
	print "4.ppid: $fields[3]\n";
	print "5.process group ID: $fields[4]\n";
	print "6.session ID: $fields[5]\n";
	print "7.number of controlling terminal: $fields[6]\n";
	print "8.process group ID in the foreground: $fields[7]\n";
	print "9.task flags: $fields[8]\n";
	print "10.number of minor page faults: $fields[9]\n";
	print "11.number of minor page faults, including those from child processes: $fields[10]\n";
	print "12.number of major page faults: $fields[11]\n";
	print "13.number of major page faults, including those from child processes: $fields[12]\n";
	print "14.CPU time in user mode (jiffies): $fields[13]\n";
	print "15.CPU time in kernel mode (jiffies): $fields[14]\n";
	print "16.CPU time in user mode (jiffies), including time from children: $fields[15]\n";
	print "17.CPU time in kernel mode (jiffies), including time from children: $fields[16]\n";
	print "18.priority: $fields[17]\n";
	print "19.niceness: $fields[18]\n";
	print "20.number of threads: $fields[19]\n";
	print "21.obsolete field: $fields[20]\n";
	print "22.time when the process started: $fields[21]\n";
	print "23.size of virtual memory space: $fields[22]\n";
	print "24.size of resident set: $fields[23]\n";
	print "25.limit of resident set size: $fields[24]\n";
	print "26.start address of the program code: $fields[25]\n";
	print "27.end address of the program code: $fields[26]\n";
	print "28.start address of the stack: $fields[27]\n";
	print "29.current value of the stack pointer: $fields[28]\n";
	print "30.current value of the instruction pointer: $fields[29]\n";
	print "31.bitmask of pending signals: $fields[30]\n";
	print "32.bitmask of blocked signals: $fields[31]\n";
	print "33.bitmask of ignored signals: $fields[32]\n";
	print "34.bitmask of signals for which a handler is set: $fields[33]\n";
	print "35.address of a kernel function where the process currently waits: $fields[34]\n";
	print "36.placeholder: $fields[35]\n";
	print "37.placeholder: $fields[36]\n";
	print "38.signal number sent to the parent when it terminates: $fields[37]\n";
	print "39.current CPU core ID: $fields[38]\n";
	print "40.realtime priority: $fields[39]\n";
	print "41.scheduling policy: $fields[40]\n";
	print "42.time the process has spent waiting for I/O: $fields[41]\n";
	print "43.guest time: $fields[42]\n";
	print "44.guest time of the children: $fields[43]\n";
}
