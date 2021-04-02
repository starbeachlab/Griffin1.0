#!/usr/bin/perl -w

use strict;

#####   ADJUST MARKED LINE BELOW IF NOT USING NAMD FORMATTED INPUT FILES !!!  #####


if( scalar( @ARGV) < 5){
    print "\nUSAGE:\n\n";
    print "extract_lines.pl  COORDINATES_INPUT_FILE  FORCE_INPUT_FILE  OUTPUT  n-blocks:[FIRST_LINE LAST_LINE] (inclusive, starting from line number 1!)\n\n";
    die "\nthat's it, bye\n\n";
}

open( INA, $ARGV[0]) or die "input file $ARGV[0] not opened";

open( INB, $ARGV[1])  or die "input file $ARGV[1] not opened";

open( OUT, ">>".$ARGV[2]) or die "output file $ARGV[2] not opened";

my $line;

my $i = 1;
my $w = 0;
for( my $c = 0; $c < 0.5 * (scalar( @ARGV)-3);$c++){

    while(  $i < $ARGV[3+2*$c]){
	$line = <INA>;
	$line = <INB>;
	$i++;
    }
    
    while( $i <= $ARGV[4+2*$c]){
	
	$line = <INA>;
	$line = trim( $line);
	my @coo = split /\s+/,$line;
	
	$line = <INB>;
	$line = trim( $line);
	my @forces = split /\s+/,$line;
	
	#########################################################################
	######    ADJUST HERE IF NOT USING NAMD FORMATTED INPUT FILES    ########
	print OUT "graphics $w cone {$coo[2] $coo[3] $coo[4]} {".($coo[2]+$forces[2])." ".($coo[3]+$forces[3])." ".($coo[4]+$forces[4])."} radius 0.12\n";
	#########################       END       ###############################
	#########################################################################

	$w++;
	$i++;
    }
}


close INA;
close INB;
close OUT;


sub trim {
    my @out = @_;
    for (@out) {
        s/^\s+//;          # trim left
        s/\s+$//;          # trim right
    }
    return @out == 1
              ? $out[0]   # only one to return
              : @out;     # or many
}
