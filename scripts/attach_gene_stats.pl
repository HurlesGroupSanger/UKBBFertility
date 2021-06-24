#!/usr/bin/env perl

use strict;
use warnings;
use Math::BigFloat ':constant';
use List::Util qw(product);

open (PLI, "<rawdata/genelists/gnomad.v2.1.1.lof_metrics.by_gene.txt") || die "Cannot open file: $!";
open (SHET, "<rawdata/genelists/shet.hgnc.txt") || die "Cannot open file: $!";
open (SHETOLD, "<rawdata/genelists/shet.cassa.hgnc.txt") || die "Cannot open file: $!";
open (GENELISTS, "<rawdata/genelists/gene_lists.txt") || die "Cannot open file: !";

## Get tissue-specific lists
opendir my $dir, "rawdata/genelists/tissues/" or die "Cannot open directory: $!";
my @lists = <GENELISTS>;
close GENELISTS;
my @files = grep { /High/ } readdir $dir;
closedir $dir;

## Make sure all file lists are here:
push (@lists, @files);

my %pli;

foreach (<PLI>) {

	chomp $_;
	my @data = split("\t", $_);
	if ($data[20] ne 'NA') {
		$pli{$data[63]} = $data[20];
	}
	
}

close PLI;

my %shet;

foreach (<SHET>) {

	chomp $_;
	my @data = split("\t", $_);
	$shet{$data[0]} = $data[2];
	
}

close SHET;

my %shet_old;

foreach (<SHETOLD>) {

	chomp $_;
	my @data = split("\t", $_);
	$shet_old{$data[0]} = $data[2];
	
}

close SHETOLD;

my %geneLists;

foreach (@lists) {

	chomp $_;
	
	if ($_ =~ /(\S+)\.txt/) {
		my $id = $1;
		if ($_ =~ /High/) {
			open (LIST, "<tissues/$_") || die "Cannot open file: $!";
		} else {
			open (LIST, "<$_") || die "Cannot open file: $!";
		}
		foreach my $gene (<LIST>) {

			chomp $gene;
			$geneLists{$id}{$gene} = 1;
			
		}

		close LIST;
		
	}

}

close GENELISTS;

my %dels;
my %dups;

open (ATTACH, "<$ARGV[0]") || die "Cannot open file: $!";

foreach my $line (<ATTACH>) {

	chomp $line;

	next if ($line =~ /^\#/);
	
	my @data = split("\t", $line);
	
	my $loc = $data[1];
	my $prim = $data[17];
	my $csq = $data[6];
	my $ct = $data[2];
	my $gene = $data[3];
	
	my $chr;
	my $start;
	my $end;
	if ($loc =~ /(\d+)\:(\d+)\-(\d+)/) {
		$chr = $1;
		$start = $2;
		$end = $3;
	} else {
		die "Couldn't parse location: $!";
	}

	if ($ct eq 'deletion' & !exists $dels{$loc}) {
		$dels{$loc} = {plis => [], shets => [], oldshets => [], genes => [], chr => $chr, start => $start, end => $end};
	} elsif ($ct eq 'duplication' & !exists $dups{$loc}) {
		$dups{$loc} = {plis => [], shets => [], oldshets => [], genes => [], chr => $chr, start => $start, end => $end};
	}

	my $pliScore = -1;
	if (exists $pli{$gene} & $prim eq 'YES') {
		$pliScore = $pli{$gene};
	}
	my $shetScore = -1;
	if (exists $shet{$gene} & $prim eq 'YES') {
		$shetScore = $shet{$gene};
	}
	my $oldshetScore = -1;
	if (exists $shet_old{$gene} & $prim eq 'YES') {
		$oldshetScore = $shet_old{$gene};
	}
	
	if ($pliScore != -1 || $shetScore != -1) {
		if ($csq =~ /transcript_ablation/ || ($csq =~ /coding_sequence_variant/ && $csq =~ /feature_truncation/)) {
			push (@{$dels{$loc}{genes}}, $gene);
			push (@{$dels{$loc}{plis}}, $pliScore);
			push (@{$dels{$loc}{shets}}, $shetScore);
			push (@{$dels{$loc}{oldshets}}, $oldshetScore);
		} elsif ($csq =~ /transcript_amplification/ || ($csq =~ /coding_sequence_variant/ && $csq =~ /feature_elongation/)) {
			push (@{$dups{$loc}{genes}}, $gene);
			push (@{$dups{$loc}{plis}}, $pliScore);
			push (@{$dups{$loc}{shets}}, $shetScore);
			push (@{$dups{$loc}{oldshets}}, $oldshetScore);
		}
	}
}


my @header = ('chr', 'start', 'end', 'ct', 'genes', 'plis', 'shets', 'highPLI', 'highsHET', 'product_sHET','product_sHET_old');
foreach my $list (sort keys %geneLists) {
	push (@header, "product_sHET_no_" . $list);
}

if ($ARGV[0] =~ /001/) {
	print join("\t", @header) . "\n";
}

printResults(\%dels, "DEL");
printResults(\%dups, "DUP");

sub printResults {

	my ($hash, $type) = @_;
	
	foreach my $sites (sort keys %$hash) {
	
		my $genes = join(",", @{$$hash{$sites}{genes}});
		my $pliScores = join(",", @{$$hash{$sites}{plis}});
		my $shets = join(",", @{$$hash{$sites}{shets}});

		## Calculate product shet:
		my @printer = ($$hash{$sites}{chr},$$hash{$sites}{start},$$hash{$sites}{end},$type);
		if ($pliScores) {
			push (@printer, ($genes,$pliScores, $shets));

			my $highpLI = highGenes(\@{$$hash{$sites}{plis}}, 0.9);
			push (@printer, $highpLI);
			
			my $highsHET = highGenes(\@{$$hash{$sites}{shets}}, 0.15);
			push (@printer, $highsHET);
			
			my $prod = calcShet(\@{$$hash{$sites}{genes}},\@{$$hash{$sites}{shets}});
			push (@printer, $prod);

			my $oldProd = calcShet(\@{$$hash{$sites}{genes}},\@{$$hash{$sites}{oldshets}});
			push(@printer, $oldProd);
			
			foreach my $list (sort keys %geneLists) {
				$prod = calcShetList(\@{$$hash{$sites}{genes}},\@{$$hash{$sites}{shets}},\%{$geneLists{$list}});
				push (@printer, $prod);
			}
		} else {
			push (@printer, (0, 0, "NA", "NA", "NA"));
			
			push (@printer, (0, 0)); # product_sHET
			
			foreach my $list (sort keys %geneLists) {
				push (@printer, 0);
			}
		}

		print join("\t", @printer) . "\n";
		
	}

}

sub calcShet {

	my ($genes, $shets) = @_;

	my @g = @$genes;
	my @s = @$shets;

	my @to_combine;
	
	for (my $x = 0; $x < scalar(@g); $x++) {
		
		if ($s[$x] != -1) {

			push (@to_combine, 1 - $s[$x]);

		}
		
	}
	
	if (scalar(@to_combine) == 0) {
		return(0);
	} else {
		return(1 - product(@to_combine));
	}
	
}

sub calcShetList {

	my ($genes, $shets, $exclude) = @_;

	my @g = @$genes;
	my @s = @$shets;
	my %e = %$exclude;

	my @to_combine;
	
	for (my $x = 0; $x < scalar(@g); $x++) {
		
		if (!exists $e{$g[$x]} && $s[$x] != -1) {

			push (@to_combine, 1 - $s[$x]);

		}
		
	}
	
	if (scalar(@to_combine) == 0) {
		return(0);
	} else {
		return(1 - product(@to_combine));
	}

}

sub highGenes {

	my ($list, $thresh) = @_;

	my @l = @$list;

	my $high = 0;

	foreach my $score (@l) {

		if ($score >= $thresh) {
			$high++;
		}
	}
		
	return ($high);	

}
