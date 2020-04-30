#!/usr/bin/env perl

use strict;
use warnings;

open (MHQ, "<rawdata/phenofiles/mhq.data.txt") || die "Cannot open file: $!";
open (OUT, ">rawdata/phenofiles/mhq.data.processed.txt") || die "Cannot make file: $!";

my %results;

my $eid;

foreach my $indv (<MHQ>) {

	chomp $indv;
	my @data = split("\t", $indv);

	$eid = $data[0];

	if (defined $eid) {
	
		## Initialize fields I care about:
		$results{$eid} = {answered_mhq => 'NA',
						  substance => 'NA',
						  soc_anxiety => 'NA',
						  scizo => 'NA',
						  phobia => 'NA',
						  panic_attacks => 'NA',
						  ocd => 'NA',
						  bipolar => 'NA',
						  depression => 'NA',
						  eating_disorders => 'NA',
						  asd => 'NA',
						  gen_anxiety => 'NA',
						  add => 'NA',
						  agoraphobia => 'NA',
						  gen_personality => 'NA'
		};
	
		## Check if the individual responded to the MHQ in the first place
		if ($data[1] eq '0' || $data[1] eq '1') {
			## Process additions
			$results{$eid}{answered_mhq} = 1;
			my %cds;
			for (my $x = 5; $x < 21; $x++) {
				if (defined $data[$x]) {
					$cds{$data[$x]} = 1;
				}
			}
			## handle missingness codes (which are two separate question panels...);
			## also annoyingly, "none of the above" is not coded, so I am assuming if they did answer the mhq, and have a blank category here, they had none of the above...
			if (!exists $cds{'-818'}) {
				$results{$eid}{depression} = 0;
				$results{$eid}{bipolar} = 0;
				$results{$eid}{gen_anxiety} = 0;
				$results{$eid}{soc_anxiety} = 0;
				$results{$eid}{agoraphobia} = 0;
				$results{$eid}{phobia} = 0;
				$results{$eid}{panic_attacks} = 0;
				$results{$eid}{ocd} = 0;
				if (exists $cds{11}) {
					$results{$eid}{depression} = 1;
				}
				if (exists $cds{10}) {
					$results{$eid}{bipolar} = 1;
				}
				if (exists $cds{15}) {
					$results{$eid}{gen_anxiety} = 1;
				}
				if (exists $cds{1}) {
					$results{$eid}{soc_anxiety} = 1;
				}
				if (exists $cds{17}) {
					$results{$eid}{agoraphobia} = 1;
				}
				if (exists $cds{5}) {
					$results{$eid}{phobia} = 1;
				}
				if (exists $cds{6}) {
					$results{$eid}{panic_attacks} = 1;
				}
				if (exists $cds{7}) {
					$results{$eid}{ocd} = 1;
				}
			}
			
			if (!exists $cds{'-819'}) {
				$results{$eid}{eating_disorders} = 0;
				$results{$eid}{scizo} = 0;
				$results{$eid}{gen_personality} = 0;
				$results{$eid}{asd} = 0;
				$results{$eid}{add} = 0;
				if (exists $cds{12} || exists $cds{13} || exists $cds{16}) {
					$results{$eid}{eating_disorders} = 1;
				}
				if (exists $cds{2}) {
					$results{$eid}{scizo} = 1;
				}
				if (exists $cds{4}) {
					$results{$eid}{gen_personality} = 1;
				}
				if (exists $cds{14}) {
					$results{$eid}{asd} = 1;
				}
				if (exists $cds{18}) {
					$results{$eid}{add} = 1;
				}
			}
			
			if ($data[1] eq '0') {
				$results{$eid}{substance} = 0;
			} elsif ($data[1] eq '1') {
				if ($data[2] eq '1' || $data[3] eq '1' || $data[4] eq '1') {
					$results{$eid}{substance} = 1;
				} else {
					$results{$eid}{substance} = 0;
				}
			}
		}
	}
}

close MHQ;

my @header = ('eid', sort keys %{$results{$eid}});

print OUT join("\t", @header) . "\n";

foreach my $indv (keys %results) {

	my @printer = ($indv);

	foreach my $key (sort keys %{$results{$indv}}) {
		push (@printer,$results{$indv}{$key});
	}

	print OUT join("\t", @printer) . "\n";

}
