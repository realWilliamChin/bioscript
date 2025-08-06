#!/usr/bin/perl
#处理kaas注释结果文件    *.keg文件
use strict;
use warnings;
use Getopt::Long;
my %opts;
GetOptions(\%opts,'i:s','o:s');
open IN, "$opts{'i'}" || die;
open OUT, ">$opts{o}" || die;
my $A_id;
my $A_name;
my $B_id;
my $B_name;
my $C_id;
my $C_name;
while (<IN>) {
	chomp;
	if (/^(A\d+)\s(.+)$/) {
		$A_id = $1;
		$A_name = $2;
	}
	elsif (/^B  (\d+)\s(.+)$/) {
		$B_id = $1;
		$B_name = $2;
	}
#C    00960 Tropane, piperidine and pyridine alkaloid biosynthesis [PATH:ko00960]
	elsif (/^C    \d+\s(.+)\s\[PATH:(ko\d+)\]/) {
		$C_id = $2;
		$C_name = $1;
	}
#C    00535 Proteoglycans [BR:ko00535]
	elsif (/^C    \d+\s(.+)\s\[BR:(ko\d+)\]/) {
		$C_id = $2;
		$C_name = $1;
	}
#C    99974 Translation
	elsif (/^C    (\d+)\s(.+)$/) {
		$C_id = $1;
		$C_name = $2;
	}
#D      transcript_6198; K15535  PWD; phosphoglucan, water dikinase [EC:2.7.9.5]
	elsif (/^D      (.+)$/) {
		my @tmp = split /;/,$1;
		$tmp[1] =~ /(K\d+)  (.+)$/;
		my $k_id_gene = $1;
		my $k_name_gene = $2;
		print OUT $tmp[0],"\t",$C_id,":",$C_name,"\t",$B_id,":",$B_name,"\t",$A_id,":",$A_name,"\t",$k_id_gene,"\t",$k_name_gene,"\t",$tmp[-1],"\n";
	}
	else{};
}

