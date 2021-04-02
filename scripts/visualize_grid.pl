#!/usr/bin/perl -w

use strict;

#################################################
##########  ADJUST SELECTIONS  ##################

### ADJUST SPATIAL LIMITS

my $xmin = -50;
my $xmax =  50;

my $ymin = -50;
my $ymax =  50;

my $zmin = 0;
my $zmax = 0.5;


### ADJUST SELECTION OF PLOTTED GRIDPOINT TYPES
### COMMENT OUT IF NOT TO PLOT

# plot sforces (if geometric objects are used this will plot the sforces unaffected by them)
my $plot_const = ""; 
$plot_const = "true"; 

# plot interaction forces  (if geometric objects are used this will plot the interaction forces unaffected by them)
my $plot_interaction = "";
#$plot_interaction = "true";

# set this to "true" if geometric objects were used when calculating the grid !!
my $objects_are_defined = "";
$objects_are_defined = "true";


### ADJUST ONLY IF NUMBER OF CONSIDERED INTERACTION INCREASES 
### OR IF GEOMETRIC OBJECTS ARE DEFINED FOR MORE THAN ONE MOLECULE TYPE 

my $max_inter = 3;      # NUMBER OF INTERACTIONS    (coulomb+attractive_vdw+repulsive_vdw)
my $max_typem = 2;      # (NUMBER OF MOLECULE TYPES) + 1 
my $max_recurs = 4;     # (NUMBER OF INTERACTIONS) + 1       



#############  ADJUST END  ######################
#################################################

if( scalar( @ARGV) != 2){
    print "\nUSAGE:\n\n";
    print "visualize_grid.pl  INPUT_GRID  OUTPUT_VMD \n\n";
    die "bye\n";
}


open( IN, $ARGV[0]) or die "grid file not opened";
open( OUT, ">".$ARGV[1]) or die "output file not opened";


my $c = 0;
my $c_const;
my $c_inter;
my $c_recurs;
my $c_typem;

my @colors;

push( @colors, "white");
push( @colors, "green");
push( @colors, "blue");
push( @colors, "red");
push( @colors, "orange");
push( @colors, "magenta");
push( @colors, "brown");
push( @colors, "yellow");
push( @colors, "grey");
push( @colors, "tan");
push( @colors, "cyan");
push( @colors, "pink");
push( @colors, "orange2");
push( @colors, "purple");


if( $plot_const eq "true"){
    $c_const = $c;
    print OUT "mol new\n";
    print OUT "graphics $c color $colors[$c]\n";
    $c++;
}

if( $plot_interaction eq "true"){
    $c_inter = $c;

    for( my $i = $c; $i < $c + $max_inter; $i++){
	print OUT "mol new\n";
	print OUT "graphics $i color $colors[$i]\n";
    }
    $c += $max_inter;
}

if( $objects_are_defined eq "true"){

    $c_typem = $c;

    for( my $i = $c; $i < $c + $max_typem; $i++){
	print OUT "mol new\n";
	print OUT "graphics $i color $colors[$i]\n";
    }

    $c += $max_typem;

    $c_recurs = $c;

    for( my $i = $c; $i < $c + $max_recurs; $i++){
	print OUT "mol new\n";
	print OUT "graphics $i color $colors[$i]\n";
    }
    $c += $max_recurs;
}



my $line;

for( my $i = 0; $i < 6; $i++){
    $line = <IN>;
}


my @bins = split /\s+/,$line;

my $total_bins = $bins[0] * $bins[1] * $bins[2];

print "bins: $bins[0]  $bins[1]  $bins[2]  $total_bins \n";

for( my $i = 0; $i < 4; $i++){
     $line = <IN>;
}

$line = trim($line);


my @min = split /\s+/,$line;

print "min: $min[0]  $min[1]  $min[2] \n";

for( my $i = 0; $i < 4; $i++){
     $line = <IN>;
}

$line = trim($line);


my @max   = split /\s+/,$line;

print "max: $max[0]  $max[1]  $max[2] \n";
for( my $i = 0; $i < 4; $i++){
     $line = <IN>;
}

$line = trim($line);


my @delta = split /\s+/,$line;

print "delta: $delta[0]  $delta[1]  $delta[2] \n";

my $total = 0;
my @count;

my $i = 0;
my $j = 0;
my $k = 0;

my $x = $min[0];
my $y = $min[1];
my $z = $min[2];

my $p = 0;

while( <IN>){

    if( $_ =~ "GridPoint"){

	# do stuff here

	if( $_ =~ "store::Constant" && $plot_const eq "true"){
	    
	    if( $x > $xmin && $x < $xmax && $y > $ymin && $y < $ymax && $z > $zmin && $z < $zmax){
		
		<IN>;
		$line = <IN>;
		
		$p++;
		$line = trim( $line);
		my @vec = split /\s+/, $line;
		
		print OUT "graphics $c_const cone {$x $y $z} {".($x+$vec[0])." ".($y+$vec[1])." ".( $z+$vec[2])."} radius 0.12\n";
	    }
	}
	elsif( $_ =~ "store::Recursive"){

	    $line = <IN>;
	    my @cols = split /\s+/, $line;
	    my $nr = $cols[0];
	    
	    for( my $i = 0; $i < $nr; $i++){
		
		$line = <IN>;
#		    print "L: ". $line;
		
		if( $line =~ "Constant"){
		    
		    $line = <IN>;
		    $line = <IN>;
		    
		    $line = trim( $line);
		    
		    if( $x > $xmin && $x < $xmax && $y > $ymin && $y < $ymax && $z > $zmin && $z < $zmax){
			my @vec = split /\s+/, $line;
			
			print OUT "graphics $c_recurs cone {$x $y $z} {".($x+$vec[0])." ".($y+$vec[1])." ".( $z+$vec[2])."} radius 0.12\n";
		    }
		    
		    for( my $j = 0; $j < 4; $j++){
			$line = <IN>;
		    }
		}
		elsif( $line =~ "Interaction"){

		    
		    my @ll = split /\s+/, $line;
		    my $nr = $ll[-1];
		    
		    for( my $j = 0; $j < $nr; $j++){
			
			$line = <IN>;
			
			if( $x > $xmin && $x < $xmax && $y > $ymin && $y < $ymax && $z > $zmin && $z < $zmax){
			    my @vec = split /\s+/, $line;
			    
			    my $m = $c_recurs;
			    
			    
			    if( $line =~ "Electrostatic"){
				$m += 1;
			    }
			    elsif( $line =~ "Attractive"){
				$m += 2;
			    }
			    elsif( $line =~ "Repulsive"){
				$m += 3;
			    }
			    else{
				die "ERROR: expected 'Electrostatic', 'Attractive' or 'Repulsive', got $line \n";
			    }
			    
			    if ( $vec[2]*$vec[2] + $vec[3]*$vec[3] + $vec[4]*$vec[4] > 0.01){
				print OUT "graphics $m cone {$x $y $z} {".($x+$vec[2])." ".($y+$vec[3])." ".( $z+$vec[4])."} radius 0.12\n";
			    }
			}
			$line = <IN>;
			$line = <IN>;
		    }
		}
		else{
		    die "ERROR: expected ConstantGridPoint or InteractionGridPoint, got: <$line>\n";
		}
	    }
	    $line = <IN>;
#		print $line;
#		die;
	}
	elsif( $_ =~ "store::TypeMapped"){

	    my @cols = split;
	    my $nr = $cols[-1];
	    for( my $i = 0; $i < $nr; $i++){
		$line = <IN>;
		$line = <IN>;
		$line = <IN>;
		
		if( $x > $xmin && $x < $xmax && $y > $ymin && $y < $ymax && $z > $zmin && $z < $zmax){
		    $line = trim( $line);
		    my @vec = split /\s+/, $line;
		    print OUT "graphics ".($c_typem+$i)." cone {$x $y $z} {".($x+$vec[0])." ".($y+$vec[1])." ".( $z+$vec[2])."} radius 0.12\n";
		}
		$line = <IN>;
		$line = <IN>;
	    }
	}
	elsif( $_ =~ "store::Interaction" && $plot_interaction eq "true"){

	    
	    if( $x > $xmin && $x < $xmax && $y > $ymin && $y < $ymax && $z > $zmin && $z < $zmax){

		my @cols = split;
		my $nr = $cols[-1];
		my $m;

		for( my $i = 0; $i < $nr; $i++){

		    $line = <IN>;
		    my @vec = split /\s+/, $line;

		    if( $line =~ /Electrostatic/){
			$m = $c_inter;
		    }
		    elsif( $line =~ /Attractive/){
			$m = $c_inter + 1;
		    }
		    elsif( $line =~ /Repulsive/){
			$m = $c_inter + 2;
		    }

		    $line = <IN>;
		    $line = <IN>;

		    if ( $vec[2]*$vec[2] + $vec[3]*$vec[3] + $vec[4]*$vec[4] > 0.01){
			print OUT "graphics $m cone {$x $y $z} {".($x+$vec[2])." ".($y+$vec[3])." ".( $z+$vec[4])."} radius 0.12\n";
		    }
		}
	    }
#		$c += 3; 
	}
	elsif( $_ =~ "Void" || $_ =~ "Surf" || $_ =~ "Interaction" || $_ =~ "Constant"){
	}
	else{
	    die "ERROR: undefined GridPoint type: $_ \n";
	}
	

    
	# end do stuff 
    
	$total++;
	$k++;
	$z += $delta[2];
	if( $k >= $bins[2]){
	    $k = 0;
	    $j++;
	    $z = $min[2];
	    $y += $delta[1];
	    if( $j >= $bins[1]){
		$j = 0;
		$i++;
		$y = $min[1];
		$x += $delta[0];

#		print "".($bins[0] - $i )." ";
	    }
	}
#	print "$i  $j  $k \n";
    }
}

print "counters: $i  $j  $k   \n";
print "coordinates: $x  $y  $z \n";
print "totals: expected: $total_bins counted: $total \n";
print "printed: $p\n";

close IN;
close OUT;



sub WriteGridPoint {

    



}






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
