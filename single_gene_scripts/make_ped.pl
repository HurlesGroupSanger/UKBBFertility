#!/usr/bin/env perl

use strict;
use warnings;

open (SPARSE, "<$ARGV[0]") || die "Cannot open file: $!";
open (SAMPLES, "<$ARGV[1]") || die "Cannot open file: $!";
open (GENES, "<$ARGV[2]") || die "Cannot open file: $!";

my %samples;

foreach (<SPARSE>) {

	chomp $_;
	my @data = split("\t", $_);

	if (exists $samples{$data[1]}) {
		
		$samples{$data[1]}{$data[0]} = $data[2];

	} else {

		$samples{$data[1]} = {$data[0] => $data[2]}

	}

}

close SPARSE;

my @samples;
my %row_sample;
foreach my $samp (<SAMPLES>) {

	chomp $samp;
	my @s = split(" ", $samp);
	push(@samples,$s[0]);
	$row_sample{$s[0]} = $samp;

}

close SAMPLES;

my @genes;
foreach my $ge (<GENES>) {
	
	chomp $ge;
	my @g = split(" ", $ge);
	push(@genes, $g[1]);
	
}

close GENES;

open (OUTPUT, ">$ARGV[3]") || die "Cannot make file: $!";

foreach my $sample (@samples) {

	my @printer = ($row_sample{$sample});

	if (exists $samples{$sample}) {
	
		foreach my $gene (@genes) {

			if (exists $samples{$sample}{$gene}) {
				
				push(@printer, "12");
				
			} else {

				push(@printer, "11");
				
			}

		} 

	} else {
		
		my @arr = ("11") x @genes;
		push(@printer, @arr);
		
	}
	
	print OUTPUT join(" ", @printer) . "\n";
	
}
