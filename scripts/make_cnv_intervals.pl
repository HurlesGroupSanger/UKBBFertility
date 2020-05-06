#!/usr/bin/env perl

use strict;
use warnings;

open (CNVS, "-|", "tail -n +2 $ARGV[0]") || die "Cannot open file: $!";

my %intervals;

foreach my $cnv (<CNVS>) {

	chomp $cnv;
	my @data = split("\t", $cnv);

	my $ct;
	if ($data[8] < 2) {
		$ct = 'DEL';
	} else {
		$ct = 'DUP'
	}

	## Column 60 is the locus...
	if (exists $intervals{$data[48]}{$ct}) {
		next;
	}
	## This is because VEP cannot annotate excessivly long CNVs. I'm removing anything greater than 50mb right now...
	if ($data[9] > 50000000) {
		next;
	}

	if ($data[48] =~ /(\d+)\:(\d+)\-(\d+)/) {

		print "$1\t$2\t$3\t$ct\n";
		$intervals{$data[48]}{$ct} = 0;

	}

}

close CNVS;
