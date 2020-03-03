#!/bin/bash

BEDPE=$1
BED=$2

paste <(/bin/grep DEL $BEDPE | awk '$1==$4' | cut -f 1,2,3,5,6 | perl -alne '$,=" ";use List::Util qw(max min);print $F[0]."\t".min($F[1], $F[2], $F[3], $F[4])."\t".max($F[1], $F[2], $F[3], $F[4])') <(/bin/grep DEL $BEDPE | awk '$1==$4' | cut -f 1-6,11,12 | tr '\t' ':') > $BED
