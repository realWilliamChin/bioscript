#!/usr/bin/perl
use strict;
# 2014.06.30    Lianfu Chen        chenllianfu@foxmail.com
# 2019.05.15    Lianfu Chen        chenllianfu@foxmail.com

my $usage = <<USAGE;
Usage:
    $0 tRNA.out tRNA.ss > tRNA.gff3

USAGE
if (@ARGV==0){die $usage}

my (%second_structure, $ss_seqID, $ss_start, $ss_seq);
open IN, $ARGV[1] or die $!;
while (<IN>) {
    if (/^(\S+)\.trna\d+ \((\d+)-(\d+)\)/) {
        $ss_seqID = $1;
        if ($3 > $2) { $ss_start = $2 } else { $ss_start = $3 }
    }
    elsif (/^Seq: (\w+)/) { $second_structure{$ss_seqID}{$ss_start}{"seq"} = $1 }
    elsif (/^Str: (\S+)/) { $second_structure{$ss_seqID}{$ss_start}{"str"} = $1 }
}
close IN;

my (%gff3, $total_tRNA_number,$total_tRNA_length, $pseudo_num);
$pseudo_num = 0;
open IN, $ARGV[0] or die $!;
while (<IN>) {
    my @annot = split /\s+/, $_;
    if ($annot[2] =~ m/\d+/) {
        $total_tRNA_number ++;
        my ($seqID,$start,$end,$intron_start,$intron_end,$type,$score,$anti_codon,$strand);

        if ($annot[3] > $annot[2]) { $start = $annot[2]; $end = $annot[3]; $strand = "+"; $intron_start = $annot[6] - 1; $intron_end = $annot[7] + 1; }
        else { $start = $annot[3]; $end = $annot[2]; $strand = "-"; $intron_start = $annot[7] - 1; $intron_end = $annot[6] + 1; }
        $seqID = $annot[0]; $type = $annot[4]; $score = $annot[8]; $anti_codon = $annot[5];

        if (m/pseudo/) {
            $type = "pseudo_$type";
            $pseudo_num ++;
        }

        my $attribute = "Name=tRNA-$type;Anti-codon=$anti_codon;Sequence=$second_structure{$seqID}{$start}{'seq'};Structure=$second_structure{$seqID}{$start}{'str'};";
        if ($intron_start == -1) {
            $total_tRNA_length += $end - $start + 1;
            $gff3{$seqID}{$start} .= "$seqID\ttRNAscan-SE\ttRNA\t$start\t$end\t$score\t$strand\t.\tID=tRNA_number_tRNA-$type;$attribute\n";
        }
        else {
            $total_tRNA_length += $intron_start - $start + 1 + $end - $intron_end + 1;
            $gff3{$seqID}{$start} .= "$seqID\ttRNAscan-SE\ttRNA\t$start\t$intron_start\t$score\t$strand\t.\tID=tRNA_number_tRNA-$type;$attribute\n";
            $gff3{$seqID}{$start} .= "$seqID\ttRNAscan-SE\ttRNA\t$intron_end\t$end\t$score\t$strand\t.\tID=tRNA_number_tRNA-$type;$attribute\n";
        }
    }
}
close IN;

my $avg = int($total_tRNA_length / $total_tRNA_number * 100 + 0.5) / 100;
print STDERR "tRNA number\tAverage length\tTotal tRNA length\n$total_tRNA_number\t$avg\t$total_tRNA_length\n";

my @seqID = keys %gff3;
my %sort_seqID;
foreach (@seqID) { $sort_seqID{$_} = $1 if /(\d+)\D*$/ }
@seqID = sort {$sort_seqID{$a} <=> $sort_seqID{$b}} @seqID;

my $number = 0;
foreach my $seqID (@seqID) {
    foreach (sort {$a <=> $b} keys %{$gff3{$seqID}}) {
        my $out = $gff3{$seqID}{$_};
        $number ++;
        my $id = '0' x (length($total_tRNA_number) - length($number)) . $number;
        $out =~ s/number/$id/g;

        my @out = split /\n/, $out;
        my @pos;
        foreach (@out) {
            @_ = split /\t/, $_;
            push @pos, $_[3];
            push @pos, $_[4];
        }
        @pos = sort {$a <=> $b} @pos;

        @_ = split /\t/, $out[0];
        my $gene_ID = $1 if $_[8] =~ m/ID=([^\s;]+)/;
        my $name = $1 if $_[8] =~ m/Name=([^\s;]+)/;
        print "$_[0]\t$_[1]\tgene\t$pos[0]\t$pos[-1]\t$_[5]\t$_[6]\t$_[7]\t$_[8]\n";
        print "$_[0]\t$_[1]\ttRNA\t$pos[0]\t$pos[-1]\t$_[5]\t$_[6]\t$_[7]\tID=$gene_ID.tRNA;Parent=$gene_ID;Name=$name;\n";

        my $num = 0;
        foreach (@out) {
            $num ++;
            @_ = split /\t/, $_;
            print "$_[0]\t$_[1]\texon\t$_[3]\t$_[4]\t\.\t$_[6]\t$_[7]\tID=$gene_ID.tRNA.exon$num;Parent=$gene_ID.tRNA;\n";
        }

        #print $out;
    }
}

print STDERR "$total_tRNA_number tRNAs in total, including $pseudo_num Pseudo tRNAs\n";
