#!/usr/bin/perl

use strict;

sub perc {
  my ($nom,$denom) = @_;
  
  sprintf("%.2f",$nom/$denom*100);
}

die "Syntax: <BAM file> <output file> [RE digest file]" unless scalar(@ARGV)>=2;

my ($bamfile,$output,$redigest) = @ARGV;

die "Cannot find $bamfile" unless -e $bamfile;

if(defined($redigest)) {
  warn "Looking at coverage over RE digest sites is not implemented yet"
}

my $totalpairs = 0;
my $totalreads = 0;

my $mappedreads = 0;

my $mappedboth = 0;
my $mappedsingle = 0;
my $mappedneither = 0;

my $intrachromosomal = 0;
my $greater20kb = 0;

open(B,"samtools view $bamfile |") or die "Cannot read BAM file $bamfile: $!";

until(eof(B)) {
  my $read1 = <B>;
  my $read2 = <B>;
  
  chomp $read1;
  chomp $read2;
  
  #WIGTC-HISEQ:2:1101:1094:9053#AGGCAG     79      *       0       0       *       *       0       0       CAGCAAGACTAGGAAGAAACTGCATCAACTAACGAGCAAAATAACCAGCTAACATCATAATGGCAGGATCAAATTCACACATANCAATATTANNNNNNAN    ccccchhhhhhhhhhhhhhhghhhhhhhhhhhhhhhhhhhhhhhhhhghhhhhhhhhhhhhghehhhhhfhhgghhhhhhhhhB[cfggghfBBBBBBBB    XM:i:0
  
  my ($r1readname,$r1flags,$r1chr,$r1pos,$r1mapq,$r1cigar,$r1matechr,$r1matepos,undef,$r1seq,$r1qstr,$r1str) = split /\t/, $read1;
  my ($r2readname,$r2flags,$r2chr,$r2pos,$r2mapq,$r2cigar,$r2matechr,$r2matepos,undef,$r2seq,$r2qstr,$r2str) = split /\t/, $read2;

  warn "$r1readname and $r2readname appear to be different paired-end reads, file may not be sorted properly" unless $r1readname eq $r2readname;
  
  $totalpairs++;
  $totalreads+=2;
  
  $mappedreads++ unless $r1chr eq '*';
  $mappedreads++ unless $r2chr eq '*';
  
  my $bothmapped = 0;

  if( $r1chr ne '*' && $r2chr ne '*' ) {
    $mappedboth++;
    $bothmapped = 1;
  } elsif( ($r1chr ne '*' && $r2chr eq '*') || ($r1chr eq '*' && $r2chr ne '*')) {
    $mappedsingle++;
  } else {
    $mappedneither++;
  }
  
  if( $bothmapped && $r1chr eq $r2chr ) {
    $intrachromosomal++;
    
    my $dist = abs($r1pos - $r2pos);
    
    $greater20kb++ if $dist>=20000;
  }
}

close(B);

open(O,">","$output") or die "Cannot write to $output: $!";

print O "Library statistics:\n";
print O "Total Number of PETs within BAM: $totalpairs ($totalreads total reads)\n";
print O "$mappedreads reads were mapped to the genome (" . perc($mappedreads,$totalreads) . "%)\n";
print O "PET mapping:\n";
print O "\tNeither end mapped: $mappedneither (" . perc($mappedneither,$totalpairs) . "%)\n";
print O "\tOnly one end mapped (singleton): $mappedsingle (" . perc($mappedsingle,$totalpairs) . "%)\n";
print O "\tMapped on both ends: $mappedboth (" . perc($mappedboth,$totalpairs) . "%)\n";
print O "\t\tInterchromosomal: " . ($mappedboth - $intrachromosomal) . " (" .perc($mappedboth - $intrachromosomal,$mappedboth) ."%)\n";
print O "\t\tIntrachromosomal: " . ($intrachromosomal) . " (" . perc($intrachromosomal,$mappedboth) . "%)\n";
print O "\t\t<20kb between PET ends: " . ($intrachromosomal - $greater20kb) . " (" . perc($intrachromosomal - $greater20kb,$mappedboth) . "%)\n";
print O "\t\t>20kb between PET ends: " . ($greater20kb) . " (" . perc($greater20kb,$mappedboth) . "%)\n";


close(O);
