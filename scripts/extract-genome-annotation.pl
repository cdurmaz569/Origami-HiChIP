#!/usr/bin/perl

use strict;

die "Syntax: <BAM file> <output file>" unless scalar(@ARGV) >= 2;

my ($bamfile,$outputfile) = @ARGV;

open(B,"samtools view -H $bamfile |") or die "Cannot extract header information from $bamfile: $!";
open(O,">","$outputfile") or die "Cannot write to file $outputfile: $!";

while(<B>) {
  chomp;
  
  next unless /\@SQ\tSN:(\w+)\tLN:(\d+)/;
  
  my $chr = $1;
  my $chrsize = $2;
  
  print O "$chr\t$chrsize\n";
}

close(B); 
close(O);

