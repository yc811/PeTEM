#!/usr/bin/perl -w
use strict;
use Data::Dumper;

my $cutoff = 4;
my $infile = $ARGV[0];
my ($id) = $infile =~ /(.+).CGmap.gz/;
my $outCG = "$id.CG.wig";
my $outCHG = "$id.CHG.wig";
my $outCHH = "$id.CHH.wig";
#open (IF, "$infile");
open (IF, "zcat $infile|");
open (OCG, ">$outCG");
open (OCHG, ">$outCHG");
open (OCHH, ">$outCHH");

print OCG "track type=wiggle_0 name=$outCG color=38,173,84 altColor=38,173,84 viewLimits=0:1\n";
print OCHG "track type=wiggle_0 name=$outCHG color=44,180,234 altColor=44,180,234 viewLimits=0:1\n";
print OCHH "track type=wiggle_0 name=$outCHH color=249,42,54 altColor=49,42,54 viewLimits=0:1\n";
my $prechr = '';
while (<IF>) {
        chomp;
        my @ar = split("\t",$_);
        next if ($ar[7] < $cutoff);
        if ($prechr ne $ar[0]) {
                $prechr = $ar[0];
                print OCG "variableStep chrom=$prechr\n";
                print OCHG "variableStep chrom=$prechr\n";
                print OCHH "variableStep chrom=$prechr\n";
        }

        if ($ar[3] eq "CG") {
                print OCG "$ar[2]\t$ar[5]\n";
        } elsif ($ar[3] eq "CHG") {
                print OCHG "$ar[2]\t$ar[5]\n";
        } elsif ($ar[3] eq "CHH") {
                print OCHH "$ar[2]\t$ar[5]\n";
        }

#       system("~/igv/IGVTools/igvtools toTDF $outfile.wig $outfile.tdf ~/igv/IGVTools/genomes/tair10.chrom.sizes");
}
