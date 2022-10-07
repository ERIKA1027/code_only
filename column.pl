#!/usr/bin/perl
#
# reform column widthin LENMAX and put '&' at end of the line
#
$LENMAX = 131;			# maximum column number for Fortran90
while ($line = <STDIN>) {
    chop $line;
    while ( length($line) > $LENMAX ) {
	$line =~ s/(^.{1,$LENMAX})/&/;
	print $1, "&\n";
    }
    print $line, "\n";
}
