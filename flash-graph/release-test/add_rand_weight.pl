#!/usr/bin/perl

while(<STDIN>) {
	/([0-9]+)\s([0-9]+)/;
	$v = rand(100);
	print "$1 $2 $v\n";
}
