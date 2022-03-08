#!/usr/bin/perl

use strict;

my %hits;

while(<STDIN>) {
  
  chomp;
  my ($c1,$s1,$e1,$c2,$s2,$e2,$readname,$score,$strand1,$strand2,$fragchr,$frags,$frage) = split /\t/;
  
  next if $c1 ne $c2; ### skip if not on the same chromosome, cannot if on the same fragment (in a normal diploid genome) 
  next unless ($strand1 eq "-" && $strand2 eq "+"); ### only find reads in +/-
  
  my $key = "$c1\t$s1\t$e2";
  my $fragment = "$fragchr\t$frags\t$frage";
  
  $hits{$key} = {} unless defined($hits{$key});
  $hits{$key}->{$fragment} = 1;
}

for my $hit (keys(%hits)) {
  next if scalar(keys(%{$hits{$hit}})) > 1;
  
  print "$hit\n";
}
