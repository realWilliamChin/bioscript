#!/usr/bin/perl
use strict;
open IN, "Trifolium_repens_Baisanye_unigene_CDS_nr_diamond_uniq.blast" || die;
my %count;
while (<IN>) {
	chomp;
	my @tmp = split /\[/,$_;
	$tmp[-1] =~ /(.+)\]$/;
	$count{$1}++;
}
close IN;
foreach (keys %count) {
	print $_,"\t",$count{$_},"\n";
}