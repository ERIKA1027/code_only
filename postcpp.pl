#! /usr/bin/perl -w
# change line after cpp

while (<>) {
    s#'//'#//#;                 # de-quate '//' -> //
    s/^\!VECT/#VECT/;		# quate !vector specifire
    s/\!.*$//;			# remove comment
    # s/^#VECT/!VECT/;		# restore !vector
    s/^#VECT/!DEC\$ VECTOR ALWAYS/;  # SSE for intel fortran
    s/;/\n/g;			# reshape column
    print unless /^\s*$/;	# remove null lines
}


