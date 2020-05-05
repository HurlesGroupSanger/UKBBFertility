#!/usr/bin/env perl

use strict;
use warnings;

open (CHUNKS, "chunks.bed") || die "Cannot open file: $!";
my $chunk = $ARGV[0];

## First search for the chunk we intended to annotate
my $line = 1;

my $chr;
my $start;
my $end;

foreach (<CHUNKS>) {

	chomp $_;
	if ($line == $chunk) {
		
		my @data = split("\t", $_);
		$chr = $data[0];
		$start = $data[1];
		$end = $data[2];

		last;
	} else {
		$line++;
	}

}

close CHUNKS;

## Set filename to use
my $fn = $chr . "_" . $start . "_" . $end;
my $cmd;

## pull variants direct from the original plink file -- this is to avoid issues with REF alleles
$cmd = "plink --bed ../ukb_efe_chr1_v1.bed --bim ../ukb_fe_exm_chrall_v1.bim --fam ../ukb44165_efe_chr1_v1_s49959.fam --recode vcf --chr $chr --from-bp $start --to-bp $end --output-chr chrMT --out $fn --real-ref-alleles --memory 10000";

runCmd($cmd, "plink");

$cmd = "bcftools view -G -O z -o $fn.sites.vcf.gz $fn.vcf";

runCmd($cmd, "bcftools view");

## Run VEP -- Make sure to set paths here:
$cmd = "perl -I /path/to/ensembl-vep/ /path/to/ensembl-vep/vep --force_overwrite --check_existing --offline --everything --no_check_variants_order --format vcf --cache --assembly GRCh38 --dir_plugins /path/to/.vep/Plugins/loftee_hg38/ --plugin LoF,loftee_path:/path/to/vep/Plugins/loftee_hg38,human_ancestor_fa:/path/to/.vep/Plugins/loftee_hg38/hg38/human_ancestor.fa.gz,conservation_file:/path/to/.vep/Plugins/loftee_hg38/hg38/loftee.sql,gerp_bigwig:/path/to/.vep/Plugins/loftee_hg38/hg38/gerp_conservation_scores.homo_sapiens.GRCh38.bw --input_file $fn.sites.vcf.gz --output_file $fn.vep.vcf --vcf";

runCmd($cmd, "VEP");

$cmd = "bgzip $fn.vep.vcf";

runCmd($cmd, "bgzip");

$cmd = "tabix -p vcf $fn.vep.vcf.gz";

runCmd($cmd, "tabix");

$cmd = "bcftools +split-vep -f \"\%CHROM\\t\%POS\\t\%REF\\t\%ALT\\t\%Gene\\t\%Feature\\t\%CANONICAL\\t\%Consequence\\t\%LoF\\t\%EXON\\t\%INTRON\\n\" -d $fn.vep.vcf.gz > $fn.vep.sites.tsv";

runCmd($cmd, "bcftools split-vep");

sub runCmd {
	
	my ($cmd, $prog) = @_;
	my $exitcd = system($cmd);

	if ($exitcd != 0) {
		die "$prog did not run correctly: $!";
	}
}
