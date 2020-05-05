#!/usr/bin/env perl

use strict;
use warnings;

open (CHUNKS, "<chunks.bed") || die "Cannot open file: $!";

foreach (<CHUNKS>) {

	chomp $_;
	
	my @data = split("\t", $_);

	my $fn = $data[0] . "_" . $data[1] . "_" . $data[2] . ".vep.sites.tsv";
	
	if (-e $fn) {

		open (SITES, "<$fn") || die "Cannot open file: $!";

		foreach my $line (<SITES>) {
			
			chomp $line;
			print "$line\n";

		}

		close SITES;

	}

}

close CHUNKS;
	
