#!/usr/bin/perl

use strict;
use warnings;
use Data::Dumper;

open (ICD, "<rawdata/phenofiles/ICD10.data.txt") || die "Cannot open file: $!";
open (OUTPUT, ">rawdata/phenofiles/ICD10.data.processed.txt") || die "Cannot make file: $!";

my %codes = (F20 => 'scizo',
			 F200 => 'scizo',
			 F201 => 'scizo',
			 F202 => 'scizo',
			 F203 => 'scizo',
			 F204 => 'scizo',
			 F205 => 'scizo',
			 F206 => 'scizo',
			 F208 => 'scizo',
			 F209 => 'scizo',
			 F231 => 'scizo',
			 F232 => 'scizo',
			 F25 => 'scizo',
			 F250 => 'scizo',
			 F251 => 'scizo',
			 F252 => 'scizo',
			 F258 => 'scizo',
			 F259 => 'scizo',
			 F840 => 'asd',
			 F841 => 'asd',
			 F845 => 'asd',
			 F849 => 'asd',
			 F30 => 'bipolar',
			 F300 => 'bipolar',
			 F301 => 'bipolar',
			 F302 => 'bipolar',
			 F308 => 'bipolar',
			 F309 => 'bipolar',
			 F31 => 'bipolar',
			 F310 => 'bipolar',
			 F311 => 'bipolar',
			 F312 => 'bipolar',
			 F313 => 'bipolar',
			 F314 => 'bipolar',
			 F315 => 'bipolar',
			 F316 => 'bipolar',
			 F317 => 'bipolar',
			 F318 => 'bipolar',
			 F319 => 'bipolar',
			 F32 => 'depression',
			 F320 => 'depression',
			 F321 => 'depression',
			 F322 => 'depression',
			 F323 => 'depression',
			 F328 => 'depression',
			 F329 => 'depression',
			 F33 => 'depression',
			 F330 => 'depression',
			 F331 => 'depression',
			 F332 => 'depression',
			 F333 => 'depression',
			 F334 => 'depression',
			 F338 => 'depression',
			 F339 => 'depression',
			 F50 => 'eating_disorders',
			 F500 => 'eating_disorders',
			 F501 => 'eating_disorders',
			 F502 => 'eating_disorders',
			 F505 => 'eating_disorders',
			 F508 => 'eating_disorders',
			 F509 => 'eating_disorders',
			 F400 => 'agoraphobia',
			 F401 => 'soc_anxiety',
			 F42 => 'ocd',
			 F420 => 'ocd',
			 F421 => 'ocd',
			 F422 => 'ocd',
			 F428 => 'ocd',
			 F429 => 'ocd',
			 F411 => 'gen_anxiety',
			 F402 => 'phobia',
			 F60 => 'gen_personality',
			 F600 => 'gen_personality',
			 F601 => 'gen_personality',
			 F602 => 'gen_personality',
			 F603 => 'gen_personality',
			 F604 => 'gen_personality',
			 F605 => 'gen_personality',
			 F606 => 'gen_personality',
			 F607 => 'gen_personality',
			 F608 => 'gen_personality',
			 F609 => 'gen_personality',
			 F61 => 'gen_personality',
			 F410 => 'panic_attacks',
			 F900 => 'add',
			 F80 => 'developmental_disorder',
			 F800 => 'developmental_disorder',
			 F801 => 'developmental_disorder',
			 F802 => 'developmental_disorder',
			 F803 => 'developmental_disorder',
			 F809 => 'developmental_disorder',
			 F81 => 'developmental_disorder',
			 F810 => 'developmental_disorder',
			 F812 => 'developmental_disorder',
			 F819 => 'developmental_disorder',
			 F82 => 'developmental_disorder',
			 F89 => 'developmental_disorder',
			 F70 => 'developmental_disorder',
			 F700 => 'developmental_disorder',
			 F701 => 'developmental_disorder',
			 F708 => 'developmental_disorder',
			 F709 => 'developmental_disorder',
			 F71 => 'developmental_disorder',
			 F711 => 'developmental_disorder',
			 F719 => 'developmental_disorder',
			 F72 => 'developmental_disorder',
			 F729 => 'developmental_disorder',
			 F78 => 'developmental_disorder',
			 F780 => 'developmental_disorder',
			 F789 => 'developmental_disorder',
			 F79 => 'developmental_disorder',
			 F790 => 'developmental_disorder',
			 F799 => 'developmental_disorder',
			 N46 => 'infertility'
	);

## For substance abuse going to automatically add since it's such a broad set of categories:
for (my $x = 10; $x < 20; $x++) {
	for (my $y = 0; $y < 5; $y++) {
		$codes{"F$x" . "$y"} = 'substance';
	}
	for (my $y = 6; $y < 10; $y++) {
		$codes{"F$x" . "$y"} = 'substance';
	}
}

my %results;

my $eid;
my %tots;

foreach my $indv (<ICD>) {

	chomp $indv;
	my @data = split("\t", $indv);
	
	$eid = $data[0];

	if (defined $eid) { 
		$results{$eid} = {scizo => 0,
						  asd => 0,
						  bipolar => 0,
						  depression => 0,
						  eating_disorders => 0,
						  substance => 0,
						  agoraphobia => 0,
						  soc_anxiety => 0,
						  ocd => 0,
						  gen_anxiety => 0,
						  phobia => 0,
						  gen_personality => 0,
						  panic_attacks => 0,
						  add => 0,
						  infertility => 0,
						  developmental_disorder => 0};
						  	
		## Process additions
		for (my $x = 1; $x < scalar(@data); $x++) {
			if (defined $data[$x]) {
				if (exists $codes{$data[$x]}) {
					$results{$eid}{$codes{$data[$x]}} = 1;
				}
			}
		}
	}
}

my @header = ('eid', sort keys %{$results{$eid}});

print OUTPUT join("\t", @header) . "\n";

foreach my $indv (keys %results) {

	my @printer = ($indv);

	foreach my $key (sort keys %{$results{$indv}}) {
		push (@printer,$results{$indv}{$key});
	}

	print OUTPUT join("\t", @printer) . "\n";

}

close ICD;
close OUTPUT;
