#!/usr/bin/perl

use strict;

my %file_size;
my %file_start;
my $trace_file = $ARGV[0];
my %file_max_size = (
	"neostore.nodestore.db", 44040024,
	"neostore.propertystore.db", 3027489446,
	"neostore.relationshipstore.db", 2277791604,
);

#print "reconstruct $trace_file\n";

open FILE, $trace_file or die $!;

while (<FILE>) {
	/([^,]+),([0-9]+),([0-9]+),/;
	my $file_name = $1;
	my $off = $2;
	my $size = $3;
	if (defined $file_max_size{$file_name}) {
		my $max_size = $file_max_size{$file_name};
		my $loc = $off + $size;
		if ($loc > $max_size) {
#			print STDERR "access $file_name at $loc > $max_size\n";
			next;
		}
	}
	if (not defined $file_size{$file_name}) {
		$file_size{$file_name} = $off + $size;
	}
	elsif ($file_size{$file_name} < $off + $size) {
		$file_size{$file_name} = $off + $size;
	}

	if (not defined $file_start{$file_name}) {
		$file_start{$file_name} = $off;
	}
	elsif ($file_start{$file_name} > $off) {
		$file_start{$file_name} = $off;
	}
}

my $accum_size = 0;
my %file_offsets;
foreach my $key (keys %file_size) {
	my $orig = $file_size{$key};
	$file_size{$key} = ($orig + 4096) & (~4095);
	print STDERR "$key has $file_size{$key} bytes, min offset: $file_start{$key}\n";
	$file_offsets{$key} = $accum_size;
	$accum_size += $file_size{$key};
#	print "$key starts from $file_offsets{$key}\n";
}

close(FILE);

open FILE, $trace_file or die $!;

my %file_prev_off;
my %file_prev_size;
while (<FILE>) {
	/([^,]+),([0-9]+),([0-9]+),([a-zA-Z]+)/;
	my $file_name = $1;
	my $off = $2;
	my $size = $3;
	if ($4 eq "seek") {
		next;
	}
	if (defined $file_max_size{$file_name}) {
		my $max_size = $file_max_size{$file_name};
		my $loc = $off + $size;
		if ($loc > $max_size) {
#			print STDERR "access $file_name at $loc > $max_size\n";
			next;
		}
	}
	if (defined $file_offsets{$file_name}) {
		$off += $file_offsets{$file_name};
		if (not defined $file_prev_off{$file_name}) {
			$file_prev_off{$file_name} = $off;
			$file_prev_size{$file_name} = $size;
			next;
		}
		my $prev_off = $file_prev_off{$file_name};
		my $prev_size = $file_prev_size{$file_name};
		if ($off + $size == $prev_off) {
			$file_prev_off{$file_name} = $off;
			$file_prev_size{$file_name} += $size;
		}
		elsif ($prev_off + $prev_size == $off) {
			$file_prev_size{$file_name} += $size;
		}
		else {
			print "0,$file_prev_off{$file_name},$file_prev_size{$file_name},$4,\n";
			$file_prev_off{$file_name} = $off;
			$file_prev_size{$file_name} = $size;
		}
#		print "$file_name offset $file_offsets{$file_name}\n";
	}
}
