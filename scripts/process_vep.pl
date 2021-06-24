#!/usr/bin/env perl

use strict;
use warnings;

my $vep_vcf = $ARGV[0];
my $processed_vep = $ARGV[1];

my $cmd = "bcftools +split-vep -f \"\%CHROM\\t\%POS\\t\%REF\\t\%ALT\\t\%ALLELE_NUM\\t\%Gene\\t\%Feature\\t\%CANONICAL\\t\%Consequence\\t\%LoF\\t\%EXON\\t\%INTRON\\n\" -d $vep_vcf";

print "$cmd\n";

open (VEP, "-|", $cmd) || die "BCFTools failed: $!";
open (OUT, ">$processed_vep") || die "Cannot make file: $!";

foreach my $vep (<VEP>) {

	chomp $vep;
	my @data = split("\t", $vep);
	$data[0] =~ s/chr//; ## Sub out the chr annotation to make my life easier...
	my $allele_num = splice(@data, 4, 1);
	my @alleles = split(",", $data[3]);
	
	## This should print the right allele...
	$data[3] = $alleles[$allele_num - 1];
	print OUT join("\t", @data) . "\n";
	
}

close OUT;
close VEP;
