#!/usr/bin/perl

use strict;

my $addresses = "";
while (<STDIN>) {
	if (/#[0-9]+\s(0x[a-z0-9]+)\s/) {
		$addresses = $addresses . $1 . " ";
	}
}
print "$addresses\n";
system("addr2line -e $ARGV[0] $addresses");
