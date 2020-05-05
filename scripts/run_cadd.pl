#!/usr/bin/env perl

use strict;
use warnings;

my $job = $ARGV[0];;

$job = sprintf("%04d", $job);

my $cmd = "/lustre/scratch115/teams/hurles/users/eg15/INTERVAL/snv_cognition/CADD-scripts/CADD.sh -g GRCh38 -v v1.5 -o to_annotate.$job.tsv.gz to_annotate.$job.vcf";

my $out = system($cmd);

if ($out != 0) {
	die "CADD failed to run: $!";
}
