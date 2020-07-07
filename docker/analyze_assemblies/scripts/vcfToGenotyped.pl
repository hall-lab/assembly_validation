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
        print "##INFO=<ID=COV,Number=.,Type=Integer,Description=\"Contig coverage at this position\">\n";
        print "##INFO=<ID=QNAME,Number=.,Type=Integer,Description=\"Query name of variant\">\n";
        print "##INFO=<ID=QSTART,Number=.,Type=Integer,Description=\"Query start pos of variant\">\n";
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
            if ($prev_info =~ /COV=1;/) {
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
            my @cov = ();
            my @qname = ();
            my @qstart = ();
            my @prev_info_array = split(";", $prev_info);
            my @info_array = split(";", $info);
            foreach my $info_entry (c(@prev_info_array, @info_array)) {
                my ($key, $value) = split("=", $info_entry);
                if ($key eq "COV") {
                    push @cov, $value;
                }
                else if ($key eq "QNAME") {
                    push @qname, $value;
                }
                else if ($key eq "QSTART") {
                    push @qstart, $value;
                }
            }
            $prev_info = join(";", "COV=".join(",", @cov), "QNAME=".join(",", @qname), "QSTART=".join(",", @qstart));
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
            if ($prev_info =~ /COV=1;/) {
                $prev_gt = "1/1";
            }
            else {
                $prev_gt = "0/1";
            }
        }
    }
}
print join("\t", $prev_pos, $prev_id, $prev_ref, $prev_alt, $prev_qual, $prev_filter, $prev_info, $prev_format, $prev_gt)."\n";
