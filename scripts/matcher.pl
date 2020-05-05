#!/usr/bin/env perl

use strict;
use warnings;
use Getopt::Long;
use Data::Dumper;

my %args;
GetOptions(\%args,
		   "file1=s",
		   "file2=s",
		   "col1=s",
		   "col2=s",
		   "r",
		   "i",
		   "w");

if ($args{file1}) {
	open (LIST1, "<$args{file1}") || die "Cannot open file1: $!";
} else {
	die "-file1 not set: $!";
}
if ($args{file1}) {
	open (LIST2, "<$args{file2}") || die "Cannot open file2: $!";
} else {
	die "-file2 not set: $!";
}
my $col1;
my $col2;
if ($args{col1}) {
	$col1 = $args{col1};
} else {
	$col1 = 0;
}
if ($args{col2}) {
	$col2 = $args{col2};
} else {
	$col2 = 0;
}
my $rev = $args{r};
my $inc = $args{i}; ## This says keep any IDs that are not found in list 1
my $ws = $args{w}; ## This means delimit by whitespace instead of tab

my %found;

my $size = 0;

foreach (<LIST1>) {

	chomp $_;
	my @data;
	if ($ws) {
		@data = split(' ', $_);
	} else {
		@data = split("\t", $_);
	}
	my @col1_split = split(",", $col1);
	my @key;
	foreach my $col (@col1_split) {
		if ($col >= scalar(@data)) {
			die "col # for file 1 is too high\n";
		}
		push(@key, $data[$col]);
	}
	my $k = join(":", @key);
	next if ($k eq "");
	$found{$k}=$_;
	if ($size == 0) {
		$size = scalar(@data);
	} else {
		if ($size != scalar(@data)) {
			print "$_\n";
			die "Column lengths in file 1 are not identical";
		}
	}
}

close LIST1;

## Builds a filler sequence of identical column length to file 1 to add in case using -i
my @filler;
for (my $x = 0; $x < $size; $x++) {
	push(@filler, "NA");
}
my $fill = join("\t", @filler);

foreach (<LIST2>) {

	chomp $_;
	my @data;
	if ($ws) {
		@data = split(' ', $_);
	} else {
		@data = split("\t", $_);
	}
	my @col2_split = split(",", $col2);
	my @key;
	foreach my $col (@col2_split) {
		if ($col >= scalar(@data)) {
			die "col # for file 2 is too high\n";
		}
		push(@key, $data[$col]);
	}
	my $k = join(":", @key);
	if ($rev) {
		if (exists $found{$k}) {
			print "$_\t$found{$k}\n";
		} elsif ($inc) {
			print "$_\t$fill\n";
		}
	} else {
		unless (exists $found{$k}) {
			print "$_\n";
		}
	}
}

close LIST2;
