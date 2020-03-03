#!/usr/bin/env genome-perl

use strict;
use warnings;
use Data::Dumper;
use List::Util qw(max);

while (my $line = <>) {
    chomp $line;
    if ($line =~ "^#") {
        print "$line\n"
    }
    else {
        my @fields = split("\t", $line);
        my $ref_span = $fields[4]-$fields[2]+2; #add two because of padding
        my $query_span = 0;
        if ($fields[8] eq '-') {
            $query_span = $fields[13]-$fields[17];
        }
        else {
            $query_span = $fields[16]-$fields[14];
        }
        $fields[6] = $fields[6].":GAP_DIFF=".abs(max(0,$ref_span)-max(0,$query_span));
        $fields[6] = $fields[6].":R_OLAP=$ref_span";
        $fields[6] = $fields[6].":QUERY_OLAP=$query_span";
        print join("\t", @fields)."\n";
    }
}
