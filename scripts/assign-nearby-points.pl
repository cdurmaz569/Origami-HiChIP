#!/usr/bin/perl

use strict;

my $basearrayreturnlengthcheck = 18; ### update in case the size of the base return array for self-ligating bins changes

my ($binfile,$w,$p,$intfile,$outputfile) = @ARGV;

die "assign-nearby-points.pl <bin file> <window size> <local size> <input file or -> [output file]" unless scalar(@ARGV) >= 4;

my %bininfo;
my @petcounts;
my %coordtree;
my @putints;
my $fh;

open(B,"<","$binfile") or die "Cannot read bin file $binfile: $!";

while(<B>) {
  chomp;
  
  my(undef,undef,undef,$idx,$resites,$chip) = split /\t/;

  
  $bininfo{$idx} = [$resites,$chip];
}

close(B);

if( $intfile eq "-" ) {
  $fh = \*STDIN;
} else {
  open($fh,"<",$intfile) or die "Cannot read $intfile: $!";
}

my $putintidx = 0;

my $lastbinidx = -1;
my $coordtreeref;
my $bininforef;

while(<$fh>) {
  chomp;
  
  my ($idx1, $idx2, $petcount, $distance) = split /\t/;
  
  
  if($lastbinidx != $idx1) {
    $bininforef = $bininfo{$idx1};
    $coordtreeref = $coordtree{$idx1};
  }
  
  my $idxd = abs($idx2 - $idx1);
  my $secondsite = $bininfo{$idx2};
  
  my $re1 = $bininforef->[0];
  my $re2 = $secondsite->[0];
  
  my $retotal = $re1 + $re2;
  
  my $peak1 = $bininforef->[1];
  my $peak2 = $secondsite->[1];

  my $peaktotal = $peak1 + $peak2;
  
  if( $idx1 == $idx2 ) { ## skip if same bin
    push @putints, [$idx1,$idx2,$idxd,$petcount,$distance,$re1,$re2,$retotal,$peak1,$peak2,$peaktotal,-1,-1,-1,-1,-1,-1,-1];
    next;
  }
  
  $putintidx++;
  
  unless(defined($coordtreeref)) {
    $coordtreeref = {};
    $coordtreeref->{$idx2} = $putintidx;
    $coordtree{$idx1} = $coordtreeref;
  } else {
    $coordtreeref->{$idx2} = $putintidx;
  }
  
  push @putints, [$idx1,$idx2,$idxd,$petcount,$distance,$re1,$re2,$retotal,$peak1,$peak2,$peaktotal];
  push @petcounts, $petcount;
}

close($fh) unless $intfile eq "-";

my $oh;

if( $intfile eq "-" ) {
  $oh = \*STDOUT;
} else {
  open($oh,">","$outputfile") or die;
}

for my $putint (@putints) {
  
  my @parr = @{$putint};
  
  my $base = join("\t",@parr);
  my $info = "";
  
  unless(scalar(@parr)>=$basearrayreturnlengthcheck) {
    my ($mi,$mj) = ($putint->[0],$putint->[1]);
  
    my $widerl = $putint->[0]-$w;
    my $widerr = $putint->[0]+$w;
    my $widecl = $putint->[1]-$w;
    my $widecr = $putint->[1]+$w;
    
    my $localrl = $putint->[0]-$p;
    my $localrr = $putint->[0]+$p;
    my $localcl = $putint->[1]-$p;
    my $localcr = $putint->[1]+$p;

    my @wide;
    my @loc;
    my @lowerleft;
    my @vert;
    my @horz;
    
    my $rself = $coordtree{$putint->[0]};
    my $selfidx = $rself->{$putint->[1]};
    
    push @loc, $selfidx;
    push @wide, $selfidx;
    
    for( my $i = $widerl; $i <= $widerr; $i++ ) {
      next if $i < 0;
      
      my $r = $coordtree{$i};
      
      for( my $j = $widecl; $j <= $widecr; $j++ ) {
        next if $i >= $j; ## keep on the upper triangle
        next if ($i == $putint->[0] && $j == $putint->[1]);
        
        next unless(defined($r->{$j}));
        
        my $ridx = $r->{$j};

        push @wide, $ridx;
        push @loc, $ridx if( $localrl <= $i && $i <= $localrr && $localcl <= $j && $j <= $localcr );
        push @lowerleft, $ridx if( $j < $localcl && $putint->[0] <= $i && $i <= $widerr );
        push @vert, $ridx if( $mj-1 <= $j && $j <= $mj+1 && (($widerl <= $i && $i < $localrl) || ($localrr < $i && $i <= $widerr)) );
        push @horz, $ridx if( $mi-1 <= $i && $i <= $mi+1 && (($widecl <= $j && $j < $localcl) || ($localcr < $j && $j <= $widecr)) );
      }
    }
    
    #$info = "\t" . join(",",@wide) . "\t" . join(",",@loc) . "\t" . join(",",@lowerleft) . "\t" . join(",",@vert) . "\t" . join(",",@horz);
    
    my $widen = 0;
    $widen += $petcounts[$_] for(@wide);
    
    my $localn = 0;
    $localn += $petcounts[$_] for(@loc);
    
    my $lln = 0;
    $lln += $petcounts[$_] for(@lowerleft);
    
    my $vertn = 0;
    $vertn += $petcounts[$_] for(@vert);
    
    my $horzn = 0;
    $horzn += $petcounts[$_] for(@horz);
    
    $info = "\t" . $petcounts[$selfidx] . "\t$widen\t$localn\t" . ($widen-$localn) . "\t$lln\t$vertn\t$horzn";
  }
  
  print $oh "$base$info\n";
}

close($oh) unless $intfile eq "-";
