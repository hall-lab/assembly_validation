#!/usr/bin/env genome-perl

use strict;
use warnings;

my $count = 0;
while (my $line = <>) {
    chomp $line;
    my @f = split("\t", $line);
    if ($f[0] eq $f[3]) {
        if ($f[1] > $f[4]) {
            my $tempStart = $f[1];
            my $tempEnd = $f[2];
            $f[1] = $f[4];
            $f[2] = $f[5];
            $f[4] = $tempStart;
            $f[5] = $tempEnd;
            $count++;
            print STDERR $f[10]."\t".$f[6]."\n";
        }
    }
    print join("\t", @f)."\n";
}
print STDERR "Rearranged $count breakpoints\n";
