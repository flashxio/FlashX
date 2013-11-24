#!/usr/bin/perl

use strict;

my $num_args = @ARGV;
if ($num_args < 2) {
	print STDERR "init_test root_conf data_file1 [data_files ...]\n";
	exit 1;
}
my $test_roots = $ARGV[0];
my @data_files;
for (my $i = 1; $i < $num_args; $i++) {
	push(@data_files, $ARGV[$i]);
}
my $num_files = @data_files;
print "There are $num_files data files\n";
open(FILE, $test_roots) or die "can't open $test_roots: $!";

my $id = 0;
while (<FILE>) {
	if (/[0-9]+:(.+)/) {
		for (my $i = 0; $i < $num_files; $i++) {
			my $data_file = $data_files[$i];
			system("mkdir -p $1/$data_file");
			system("touch $1/$data_file/$id");
		}
		$id++;
	}
}
