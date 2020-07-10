#!/usr/bin/env genome-perl

use strict;
use warnings;
use Data::Dumper;

my $prev_pos;
my $prev_id;
my $prev_ref;
my $prev_alt;
my $prev_qual;
my $prev_filter;
my $prev_info;
my $prev_format;
my $prev_gt;

while (my $line = <>) {
    chomp $line;
    if ($line =~ "^##") {
        print $line."\n";
    }
    elsif ($line =~ "^#") {
        print "##INFO=<ID=COV,Number=1,Type=Integer,Description=\"Contig coverage at this position\">\n";
        print "##FORMAT=<ID=GT,Number=1,Type=String,Description=\"Genotype\">\n";
        print $line."\n";
    }
    else {
        my @f = split("\t", $line);
        my $pos = join("\t", $f[0], $f[1]);
        my $ref = $f[3];
        my $alt = $f[4];
        my $id = $f[2];
        my $qual = $f[5];
        my $filter = $f[6];
        my $info = $f[7];
        my $format = $f[8];
        unless (defined $prev_pos) {
            $prev_pos = $pos;
            $prev_ref = $ref;
            $prev_alt = $alt;
            $prev_id = $id;
            $prev_qual = $qual;
            $prev_filter = $filter;
            $prev_info = $info;
            $prev_format = $format;
            if ($prev_info eq "COV=1") {
                $prev_gt = "1/1";
            }
            else {
                $prev_gt = "0/1";
            }
            next;
        }
        if ($pos eq $prev_pos) {
            if ($prev_ref eq $ref and $prev_alt eq $alt) {
                $prev_gt = "1/1";
                $prev_qual = $prev_qual + $qual;
            }
            else {
                $prev_gt = "1/2";
                $prev_qual = $prev_qual + $qual;
                if (length($ref) > length($prev_ref)) {
                    my @alt_list;
                    foreach my $alt_entry (split(",", $prev_alt)) {
                        push(@alt_list, $alt_entry.substr($ref, length($prev_ref)));
                    }
                    $prev_ref = $ref;
                    $prev_alt = join(",", @alt_list);
                }
                elsif (length($ref) < length($prev_ref)) {
                    $alt = $alt.substr($prev_ref, length($ref));
                }
                $prev_alt = join(",", $prev_alt, $alt);
            }
        }
        else {
            my @uniq_alts = do { my %seen; grep { !$seen{$_}++ } split(",", $prev_alt)};
            if (scalar(@uniq_alts)>64) {
                @uniq_alts = @uniq_alts[0 .. 63];
            }
            $prev_alt = join(",", @uniq_alts);
            print join("\t", $prev_pos, $prev_id, $prev_ref, $prev_alt, $prev_qual, $prev_filter, $prev_info, $prev_format, $prev_gt)."\n";
            $prev_pos = $pos;
            $prev_ref = $ref;
            $prev_alt = $alt;
            $prev_filter = $filter;
            $prev_id = $id;
            $prev_qual = $qual;
            $prev_info = $info;
            $prev_format = $format;
            if ($prev_info eq "COV=1") {
                $prev_gt = "1/1";
            }
            else {
                $prev_gt = "0/1";
            }
        }
    }
}
print join("\t", $prev_pos, $prev_id, $prev_ref, $prev_alt, $prev_qual, $prev_filter, $prev_info, $prev_format, $prev_gt)."\n";
