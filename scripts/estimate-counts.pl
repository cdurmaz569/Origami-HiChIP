#!/usr/bin/perl

### 10, 11, 12

use strict;

die "Not enough arguments" unless (scalar(@ARGV)>=1);

my ($file,$bins) = @ARGV;

my $fh;
my %bininfo;

unless($file eq "-") {
  open($fh,"<",$file) or die "Cannot read from $file: $!";
} else {
  $fh = \*STDIN;
}

if( defined($bins) && $bins ne "" ) {
  open(B,"<","$bins") or die "Cannot read bin file $bins: $!";
  
  while(<B>) {
    chomp;
    
    my ($chr,$s,$e,$idx) = split /\t/;
    
    my $key = "$chr $s $e";
    $bininfo{$key} = $idx;
  }
  
  close(B);
}

my %c;

my $line1 = <$fh>;
my $line2 = <$fh>;

while(!eof($fh)) {
	chomp $line1;
	my @l1 = split /\t/, $line1;
	chomp $line2;
	my @l2 = split /\t/, $line2;

  if(!($l1[0] eq $l2[0] && $l1[1] == $l2[1] && $l1[2] == $l2[2] &&
      $l1[3] eq $l2[3] && $l1[4] == $l2[4] && $l1[5] == $l2[5])) {
      ## somehow we got disoriented, toss a line a move forward
      
      $line1 = $line2;
      $line2 = <$fh>;
      next;
  }

	my ($c1,$p11,$p12) = @l1[10,11,12];
	my ($c2,$p21,$p22) = @l2[10,11,12];
	
	
	### exclude chrM (although it would be nice to generalize this exclusion to arbitrary chromosomes)
	unless( $c1 =~ /chrM/ || $c2 =~ /chrM/) {

	  my $out = "";


	  if( $c1 gt $c2 ) {
		  $out = "$c2:$p21:$p22:$c1:$p11:$p12";
	  }
	  elsif( $c1 eq $c2 ) {
		  if( $p11 > $p21 ) {
	       $out = "$c2:$p21:$p22:$c1:$p11:$p12";
		  } else {
	       $out = "$c1:$p11:$p12:$c2:$p21:$p22";
		  }
	  } else {
         $out = "$c1:$p11:$p12:$c2:$p21:$p22";
	  }
	 
	  $c{$out}++;
	}
	$line1 = <$fh>;
	$line2 = <$fh>;
}


if (defined($bins) && $bins ne "") {
  
  for my $k (keys(%c)) {
  	next if $k =~ /HASH/; ## weird error...

    my ($c1,$s1,$e1,$c2,$s2,$e2) = split /:/, $k;
    
    my ($key1,$key2) = ("$c1 $s1 $e1","$c2 $s2 $e2");
    
    my ($idx1,$idx2,$counts,$d) = ($bininfo{$key1},$bininfo{$key2},$c{$k},abs($s2-$s1));
    
    $d = "NA" unless $c1 eq $c2;
    
    print "$idx1\t$idx2\t$counts\t$d\n";
  }
  
} else {
  for my $k (sort keys(%c)) {
  	next if $k =~ /HASH/; ## weird error...
  	my @arr = split /:/, $k;
  	print "$_\t" for (@arr);
  	print "$c{$k}\n";
  }
}

close($fh);
