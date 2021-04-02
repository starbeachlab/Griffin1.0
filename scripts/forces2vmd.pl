#!/usr/bin/perl -w

open( IN, $ARGV[0]) || die;

	@min = ( -120,-120,0);
	@max = ( 120, 120, 80);


print "mol new\n";
print "graphics 0 color red\n";
while( <IN>){
	next if length( $_) < 10;
	my @cols = split;
	if( 
	scalar( @cols) == 8
	    && $cols[5] * $cols[5] + $cols[6] * $cols[6] + $cols[7] * $cols[7] > 0.01
#	&& ( $cols[5] != 0 || $cols[6] != 0.0 || $cols[7] != 0.0)
#	
#	&& $min[0] <= $cols[2]
#	&& $min[1] <= $cols[3]
#	&& $min[2] <= $cols[4]
#	&& $max[0] >= $cols[2]
#	&& $max[1] >= $cols[3]
#	&& $max[2] >= $cols[4]

    ){

#         print "graphics 0 cylinder {" . $cols[ 2] . " " . $cols[3] . " " . $cols[4] . "} {" . ( $cols[5]) . " " . ( $cols[6]) . " " . ( $cols[ 7]) . "} radius 0.1\n";
#        print "graphics 0 cone {" . $cols[ 2] . " " . $cols[3] . " " . $cols[4] . "} {" . ( $cols[5] + $cols[2]) . " " . ( $cols[6] + $cols[3]) . " " . ( $cols[ 7] + $cols[4]) . "} radius 0.1\n";
        printf "graphics 0 cone {%.3f %.3f %.3f} {%.3f %.3f %.3f} radius 0.05\n", $cols[ 2],  $cols[3], $cols[4], ( $cols[5] + $cols[2]), ( $cols[6] + $cols[3]), ( $cols[ 7] + $cols[4]);
    }   
}   
