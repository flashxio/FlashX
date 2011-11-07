#!/usr/bin/perl

while (<STDIN>) {
	if (/^$/) {
		if (!($prev =~ /kernel_thread_helper/ || $prev =~ /page_fault/)) {
			print "$prev\n";
		}
	}
	$prev = $_;
}
