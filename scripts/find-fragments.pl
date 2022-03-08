#!/usr/bin/perl

use strict;

die "Syntax: <RE sequence> <FASTA file> <output BED file> <fragment file>" unless scalar(@ARGV)>=3;

my ($reseq,$fasta,$bed,$fragments) = @ARGV;

$reseq = uc $reseq;

die "$reseq needs to be a DNA sequence" unless $reseq =~ /^[ACTG]+$/;
die "Cannot find FASTA file $fasta" unless -e $fasta;

open(F,"<","$fasta") or die "Cannot read $fasta: $!";
open(O,">","$bed") or die "Cannot write to $bed: $!";
open(R,">","$fragments") or die "Cannot write to $fragments: $!";

my $relen = length $reseq;
my $chr = "";
my $seq;
my $regex = qr/$reseq/i;
my %frags;

my $revseq = $reseq;
$revseq =~ tr/CGTA/GCAT/;
my $revregex = qr/$revseq/i;

while(<F>) {
  chomp;
  
  if( /^>(.+)/ ) {
    my $chrtobe = $1;
    my @farr;

    while( $seq =~ /$regex/g ) {
      my ($s,$e) = ($-[0],$+[0]);
      print O "$chr\t$s\t$e\tRE\t1000\t+\n";

      push @farr, [$s,$e];
    }
    
    
    # doesn't seem to be a strong emphasis on the reverse case in other Hi-C pipelines
    #while ($seq =~ /$revregex/g ) {
    #  my ($s,$e) = ($-[0],$+[0]);
    #
    #  print O "$chr\t$s\t$e\tRE\t1000\t-\n";
    #
    #  push @farr, [$s,$e];
    #}
    
    @farr = sort {$a->[0] <=> $b->[0]} @farr;

    $frags{$chr} = \@farr;
    $chr = $chrtobe;
    $seq = "";
    
    next;
  }
  
  $seq .= $_;
}

if($chr ne "")  { # sanity check
  my @farr;
    
  while( $seq =~ /$regex/g ) {
    my ($s,$e) = ($-[0],$+[0]);
    print O "$chr\t$s\t$e\tRE\t1000\t+\n";

    push @farr, [$s,$e];
  }
    
    
  # doesn't seem to be a strong emphasis on the reverse case in other Hi-C pipelines
  #while ($seq =~ /$revregex/g ) {
  #  my ($s,$e) = ($-[0],$+[0]);
  #
  #  print O "$chr\t$s\t$e\tRE\t1000\t-\n";
  #
  #  push @farr, [$s,$e];
  #}
  
  @farr = sort {$a->[0] <=> $b->[0]} @farr;
    
  $frags{$chr} = \@farr;
}

close(F);
close(O);

for my $c (keys(%frags)) {
  my @a = @{$frags{$c}};
  my $x = scalar(@a);
  
  for( my $i = 1; $i < $x; $i++) {
    ### BED file format is annoying in that it represents a half-open interval with the end being open
    my $b = $a[$i-1]->[1];
    my $e = $a[$i]->[0];
    
    next if ($e-$b) < $relen; ## skip if the fragment is just two RE sites back to back
    print R "$c\t$b\t$e\tFragment\t1000\t+\n";
  }
}

close(R);
