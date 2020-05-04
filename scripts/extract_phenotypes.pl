#!/usr/bin/env perl

use strict;
use warnings;

open (HTML, "<../rawdata/UKBBPhenoFiles/ukb40103.html") || die "Cannot open file: $!";
open (OUT, ">../rawdata/UKBBPhenoFiles/pheno_fields") || die "Cannot make file: $!";

foreach my $line (<HTML>) {

  if ($line =~ /field\.cgi\?id=(\d*).*rowspan=[^>]*>([^<]+)/) {
    print OUT "$1 $2\n";
  }

}

close OUT;
close HTML;

open (FIELDS, "<../rawdata/UKBBPhenoFiles/fields_to_extract.txt") || die "Cannot open file: $!";

my $cols = '$' . "1" . '"\t"';

foreach my $line (<FIELDS>) {

  chomp $line;
  next if ($line =~ /^#/);
  my $a = $line;
  my $regex = "(\\d*)</td><td>.*field\\.cgi\\?id=$a\\\".*>[1-9]+[0-9]*</td>";
  
  open (HTML, "../rawdata/UKBBPhenoFiles/ukb40103.html") || die "Cannot open file: $!";
  foreach my $html (<HTML>) {
    #if ($html =~ /(\d*)<\/td><td>.*field\.cgi\?id=21022".*>[1-9]+[0-9]*<\/td>/) {
    if ($html =~ /$regex/) {
      $cols .= '$' . ($1 + 1) . '"\t"';
    }
  }
  close HTML;

}

close FIELDS;

$cols = substr($cols, 0, length($cols) - 4);

my $cmd = "awk -F\"\\t\" \'\{print $cols\}\' ../rawdata/UKBBPhenoFiles/ukb40103.txt > ../rawdata/UKBBPhenoFiles/ukbb_phenotypes.txt";
print "$cmd\n";

system($cmd);
