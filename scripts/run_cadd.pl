#!/usr/bin/env perl

use strict;
use warnings;

## This grabs all of the VCF files from the bed file of vcf locations
open (LIST, "-|", "zcat ../ukbb_200k.locations.bed.gz") || die "Cannot read sites file: $!";

## And then picks one of them based on the current lsf job number. This will (obviously) not work if you use something other than lsf
my $job = $ENV{'LSB_JOBINDEX'};
my $curr = 1;

my $chr;
my $start;
my $file;

## Grab the Nth file according to $job and get the chr, start, stop information
foreach my $line (<LIST>) {

	chomp $line;
	if ($curr == $job) {
		my @data = split("\t", $line);
		$chr = $data[0];
		$start = $data[1];
		$file = $data[3];
		last;
	} else {
		$curr++;
	}
}

## And get the corresponding VQSR file we already annotated
close LIST;
$chr =~ s/chr//;
print "../ukbb_200k.locations.vqsr/$chr" . "_$start.tsv\n";
open (VQSR, "<../ukbb_200k.locations.vqsr/$chr" . "_$start.tsv") || die "Cannot open file: $!";
open (OUT, ">$chr" . "_$start.vcf") || die "Cannot make file: $!";

my @to_annotate;

## And using the VQSR file to generate a CADD input file
foreach (<VQSR>) {
	chomp $_;
	my @data = split("\t", $_);
	$data[0] =~ s/chr//;
	if (length($data[2]) > 1 || length($data[3]) > 1) {
		my $cadd_str = "$data[0]\t$data[1]\t.\t$data[2]\t$data[3]";
		print OUT "$cadd_str\n";
		push(@to_annotate, $cadd_str);
	}
}

close OUT;

## Then run CADD...
my $cmd = "CADD.sh -o $chr" . "_$start.tsv.gz $chr" . "_$start.vcf";
my $err = system($cmd);

if ($err != 0) {
	die "CADD failed: $!";
}

## This takes the CADD output and matches it to the proper allele from the VQSR file above.
## This is crucial as CADD does some weird stuff when left/right correcting InDels
open (CADD, "-|", "zcat $chr" . "_$start.tsv.gz") || die "Cannot open file: $!";
open (FINAL, ">$chr" . "_$start.processed.tsv") || die "Cannot make file: $!";

my @cadd = <CADD>;
close CADD;

for (my $x = 2; $x < scalar(@cadd); $x++) {

	my $line = $cadd[$x];
	chomp $line;
	my @data = split("\t", $line);
	my @og_data = split("\t",$to_annotate[$x - 2]);

	print FINAL "$og_data[0]\t$og_data[1]\t$og_data[3]\t$og_data[4]\t$data[4]\t$data[5]\n";

}

close FINAL;
	
