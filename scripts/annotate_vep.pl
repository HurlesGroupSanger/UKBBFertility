#!/usr/bin/env perl

use strict;
use warnings;

my $job = $ARGV[0];

my $file = sprintf("%03d", $job);

# Old command for local VEP
my $cmd = "perl -I /path/to/vep_dir/ /path/to/vep_dir/vep --canonical --ccds --regulatory --numbers --force_overwrite --check_existing --offline --cache --dir /path/to/vep_cache/ --fasta /path/to/hg19_reference/Homo_sapiens.GRCh37.75.dna.primary_assembly.fa --dir_plugins --format ensembl -i cnv_for_vep.$file --tab -o vep_annotation.$file.txt";

my $ret = system($cmd);

if ($ret != 0) {

	die "VEP failed: $!";

}
