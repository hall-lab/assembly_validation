#!/usr/bin/env genome-perl

use strict;
use warnings;
use Data::Dumper;

my $a = $ARGV[0];
my $b = $ARGV[1];
my $c = $ARGV[2];
my $repType = $ARGV[3];
my $aName = $ARGV[4];
my $bName = $ARGV[5];
my $cName = $ARGV[6];

open(A, "<$a");
open(B, "<$b");
open(C, "<$c");

my %a;
my %b;
my %c;
my %aAndb;
my %bAndc;
my %aAndc;
my %aAndbAndc;

while (my $line = <A>) {
    chomp $line;
    my @fields = split("\t", $line);
    my $key = join(":", @fields[-3..-1]);
    $a{$key} = 1;
}
while (my $line = <B>) {
    chomp $line;
    my @fields = split("\t", $line);
    my $key = join(":", @fields[-3..-1]);
    $b{$key} = 1;
    if (defined $a{$key} and $a{$key} == 1) {
        $aAndb{$key} = 1;
    }
}
while (my $line = <C>) {
    chomp $line;
    my @fields = split("\t", $line);
    my $key = join(":", @fields[-3..-1]);
    $c{$key} = 1;
    if (defined $a{$key} and $a{$key} == 1) {
        $aAndc{$key} = 1;
    }
    if (defined $b{$key} and $b{$key} == 1) {
        $bAndc{$key} = 1;
    }
    if (defined $aAndb{$key} and $aAndb{$key} == 1) {
        $aAndbAndc{$key} = 1;
    }
}

close A;
close B;
close C;

print keys(%a)."\n";
print keys(%b)."\n";
print keys(%c)."\n";
print keys(%aAndb)."\n";
print keys(%bAndc)."\n";
print keys(%aAndc)."\n";
print keys(%aAndbAndc)."\n";
