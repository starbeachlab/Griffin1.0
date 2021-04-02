#!/usr/bin/perl                                                                

my $intervall = 10;

if( scalar( @ARGV) == 0){
    print "\nUSAGE:\n";
    print "monitor.pl PID\n";
    print "\ninternal parameter: \$intervall: sleeping time\n\n";
    die;
}

$control = 1;

$id = $ARGV[0];


while( $control){

    system( "pmap $id | grep otal > tmp_monitor.dat");

    open( IN, "tmp_monitor.dat");
    @lines = <IN>;
    close( IN);
    @cols = split /\s+/, $lines[0];
    if( length( $cols[2]) == 0 ){
    	$control = 0;
    }
    else{
    	print substr( $cols[2], 0, length( $cols[2]) -1)."\n";
    }
    sleep $intervall;
}

#system( "rm tmp_monitor.dat");
