#!/usr/bin/perl

# This script is used to pack the library in a tar.gz file for redistribution.

use strict;

if (@ARGV < 2) {
	print STDERR "pack pack_name version\n";
	exit;
}

my $pack_name = $ARGV[0];
my $version = $ARGV[1];

my @dirs = (
	'conf',
	'COPYING',
	'include',
	'libsafs',
	'README',
	'test',
	'utils',
);

my $rand = rand();
my $tmp_dir = "/tmp/$pack_name-v$version.$rand";
my $pack_dir = "/$tmp_dir/$pack_name-v$version";

system("mkdir -p $pack_dir");
for (my $i = 0; $i < @dirs; $i++) {
	system("cp -R $dirs[$i] $pack_dir/$dirs[$i]");
}

system("tar -zcvf $pack_name-v$version.tar.gz -C $tmp_dir $pack_name-v$version");
