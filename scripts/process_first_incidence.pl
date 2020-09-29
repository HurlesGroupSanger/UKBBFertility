#!/usr/bin/env perl

use strict;
use warnings;

open (CODES, "rawdata/phenofiles/pheno_fields") || die "Cannot open file: $!";

my %cols;
my $col = 1;

foreach (<CODES>) {

	chomp $_;
	if ($_ =~ /(\d+) ([\S\s]+)/) {
		my $code = $1;
		my $desc = $2;

		if ($desc =~ /Date ([A-Z]\d{2})/) {
			$desc = $1 . ".date";
			$cols{$col} = $desc;
			$col++;
		} elsif ($desc =~ /Source[\S\s]+([A-Z]\d{2})/) {
			$desc = $1 . ".source";
			$cols{$col} = $desc;
			$col++;
		} elsif ($desc eq "GP clinical event records") {
			$desc = "num.recs";
			$cols{$col} = $desc;
			$col++;
		}
			
	}

}

close CODES;

open (PROC, ">rawdata/phenofiles/processed_fi.txt") || die "Cannot make file: $!";
open (VALID, ">rawdata/phenofiles/valid_fi_indvs.txt") || die "Cannot make file: $!";

open (FIRST, "<rawdata/phenofiles/fi_phenotypes.txt") || die "Cannot open file: $!";

print PROC "eid\tcode\tdate\tsource\n";
print VALID "eid\tnum.gp.codes\n";

foreach (<FIRST>) {
	
	next if ($_ =~ /eid/);
	my @data = split("\t", $_);
	chomp $data[scalar(@data) - 1];
	
	my $eid = $data[0];
	my $num_reqs = $data[1];
	print VALID "$eid\t$num_reqs\n";
	if ($num_reqs) {
		for (my $x = 2; $x < scalar(@data); $x+=2) {
			my $source = $x + 1;
			if ($data[$x] ne "") {
				if ($cols{$x} =~ /([A-Z]\d{2})/) {
					print PROC "$eid\t$1\t$data[$x]\t$data[$source]\n";
				} else {
					print die "Code $eid - $cols{$x}\n$_\n";
				}
			}
		}
	}
}

close FIRST;
