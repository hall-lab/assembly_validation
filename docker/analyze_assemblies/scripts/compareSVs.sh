#!/bin/bash
set -exo pipefail

REF_SV_VCF=$1
REF_SV_BED=$2
BEDPE=$3
OUTPUT_PREFIX=$4
SAMPLE=$5
REFNAME=$6
SEGDUP_BED=$7
STR_BED=$8

SVTOOLS=/opt/hall-lab/python-2.7.15/bin/svtools
PERL=/usr/bin/perl
GREP=/bin/grep
BEDTOOLS=/opt/hall-lab/bedtools
BEDPETOBED=/opt/hall-lab/scripts/bedpetobed.sh
COUNT_REPETITIVE=/opt/hall-lab/scripts/countRepetitive.pl

REF_SV_BEDPE=$OUTPUT_PREFIX.$REFNAME.padded.bedpe
$SVTOOLS vcftobedpe -i $REF_SV_VCF -o $REF_SV_BEDPE.tmp
cat <($GREP "^#" $REF_SV_BEDPE.tmp) <( $GREP -v "^#" $REF_SV_BEDPE.tmp | $PERL -ape '$F[1] -= 1; $F[2]+=1; $F[4] -= 1; $F[5] += 1; $_ = join("\t", @F)."\n"') | awk 'BEGIN{FS="\t"; OFS=FS}; {if ($11=="INS") {$5=$2; $6=$3} print $0}'  > $REF_SV_BEDPE
$SVTOOLS bedpetobed12 -i <($GREP -w PASS $REF_SV_BEDPE) -o $OUTPUT_PREFIX.$REFNAME.forIGV.bed -n $SAMPLE.$REFNAME
$BEDTOOLS pairtobed -a $REF_SV_BEDPE -b $STR_BED -type either > $OUTPUT_PREFIX.$REFNAME.str.bedpe
$BEDTOOLS pairtobed -a <($GREP -w PASS $REF_SV_BEDPE) -b $STR_BED -type either > $OUTPUT_PREFIX.$REFNAME.str.PASS.bedpe
$BEDTOOLS pairtobed -a $REF_SV_BEDPE -b $SEGDUP_BED -type either > $OUTPUT_PREFIX.$REFNAME.allSegDup.bedpe
$BEDTOOLS pairtobed -a <($BEDTOOLS pairtobed -a $REF_SV_BEDPE -b $STR_BED -type neither) -b $SEGDUP_BED -type either > $OUTPUT_PREFIX.$REFNAME.segDup.bedpe
$BEDTOOLS pairtobed -a <($BEDTOOLS pairtobed -a <($GREP -w PASS $REF_SV_BEDPE) -b $STR_BED -type neither) -b $SEGDUP_BED -type either > $OUTPUT_PREFIX.$REFNAME.segDup.PASS.bedpe
$BEDTOOLS pairtobed -a <($BEDTOOLS pairtobed -a $REF_SV_BEDPE -b $STR_BED -type neither) -b $SEGDUP_BED -type neither > $OUTPUT_PREFIX.$REFNAME.nonRep.bedpe
$BEDTOOLS pairtobed -a <($BEDTOOLS pairtobed -a <($GREP -w PASS $REF_SV_BEDPE) -b $STR_BED -type neither) -b $SEGDUP_BED -type neither > $OUTPUT_PREFIX.$REFNAME.nonRep.PASS.bedpe

echo "total $REFNAME calls" > $OUTPUT_PREFIX.counts.txt
$GREP -v "^#" $REF_SV_BEDPE | wc -l >> $OUTPUT_PREFIX.counts.txt
echo "total $REFNAME calls PASS" >> $OUTPUT_PREFIX.counts.txt
$GREP -v "^#" $REF_SV_BEDPE | $GREP -w PASS | wc -l >> $OUTPUT_PREFIX.counts.txt
echo "total nonRep $REFNAME calls" >> $OUTPUT_PREFIX.counts.txt
$GREP -v "^#" $OUTPUT_PREFIX.$REFNAME.nonRep.bedpe | cut -f 1-6 | sort -u | wc -l >> $OUTPUT_PREFIX.counts.txt
echo "total nonRep PASS $REFNAME calls" >> $OUTPUT_PREFIX.counts.txt
$GREP -v "^#" $OUTPUT_PREFIX.$REFNAME.nonRep.PASS.bedpe | cut -f 1-6 | sort -u | wc -l >> $OUTPUT_PREFIX.counts.txt
echo "total str $REFNAME calls" >> $OUTPUT_PREFIX.counts.txt
$GREP -v "^#" $OUTPUT_PREFIX.$REFNAME.str.bedpe | cut -f 1-6 | sort -u | wc -l >> $OUTPUT_PREFIX.counts.txt
echo "total str PASS $REFNAME calls" >> $OUTPUT_PREFIX.counts.txt
$GREP -v "^#" $OUTPUT_PREFIX.$REFNAME.str.PASS.bedpe | cut -f 1-6 | sort -u | wc -l >> $OUTPUT_PREFIX.counts.txt
echo "total segDup $REFNAME calls" >> $OUTPUT_PREFIX.counts.txt
$GREP -v "^#" $OUTPUT_PREFIX.$REFNAME.segDup.bedpe | cut -f 1-6 | sort -u | wc -l >> $OUTPUT_PREFIX.counts.txt
echo "total segDup PASS $REFNAME calls" >> $OUTPUT_PREFIX.counts.txt
$GREP -v "^#" $OUTPUT_PREFIX.$REFNAME.segDup.PASS.bedpe | cut -f 1-6 | sort -u | wc -l >> $OUTPUT_PREFIX.counts.txt

$BEDTOOLS pairtopair -type both -a $OUTPUT_PREFIX.$REFNAME.nonRep.bedpe -b $BEDPE -is -slop 50 > $OUTPUT_PREFIX.compared.$REFNAME.nonRep.bedpe
echo "pairtopair for $REFNAME nonRep vs combined calls, 50 bp slop" >> $OUTPUT_PREFIX.counts.txt
cut -f 1-6 $OUTPUT_PREFIX.compared.$REFNAME.nonRep.bedpe | sort -u | wc -l >> $OUTPUT_PREFIX.counts.txt
$BEDTOOLS pairtopair -type both -a $OUTPUT_PREFIX.$REFNAME.nonRep.PASS.bedpe -b $BEDPE -is -slop 50 > $OUTPUT_PREFIX.compared.$REFNAME.nonRep.PASS.bedpe
echo "pairtopair for $REFNAME nonRep PASS vs combined calls, 50 bp slop" >> $OUTPUT_PREFIX.counts.txt
cut -f 1-6 $OUTPUT_PREFIX.compared.$REFNAME.nonRep.PASS.bedpe | sort -u | wc -l >> $OUTPUT_PREFIX.counts.txt
$BEDTOOLS pairtopair -type both -a $OUTPUT_PREFIX.$REFNAME.str.bedpe -b $BEDPE -is -slop 50 > $OUTPUT_PREFIX.compared.$REFNAME.str.bedpe
echo "pairtopair for $REFNAME str vs combined calls, 50 bp slop" >> $OUTPUT_PREFIX.counts.txt
cut -f 1-6 $OUTPUT_PREFIX.compared.$REFNAME.str.bedpe | sort -u | wc -l >> $OUTPUT_PREFIX.counts.txt
$BEDTOOLS pairtopair -type both -a $OUTPUT_PREFIX.$REFNAME.str.PASS.bedpe -b $BEDPE -is -slop 50 > $OUTPUT_PREFIX.compared.$REFNAME.str.PASS.bedpe
echo "pairtopair for $REFNAME str PASS vs combined calls, 50 bp slop" >> $OUTPUT_PREFIX.counts.txt
cut -f 1-6 $OUTPUT_PREFIX.compared.$REFNAME.str.PASS.bedpe | sort -u | wc -l >> $OUTPUT_PREFIX.counts.txt
$BEDTOOLS pairtopair -type both -a $OUTPUT_PREFIX.$REFNAME.segDup.bedpe -b $BEDPE -is -slop 50 > $OUTPUT_PREFIX.compared.$REFNAME.segDup.bedpe
echo "pairtopair for $REFNAME segDup vs combined calls, 50 bp slop" >> $OUTPUT_PREFIX.counts.txt
cut -f 1-6 $OUTPUT_PREFIX.compared.$REFNAME.segDup.bedpe | sort -u | wc -l >> $OUTPUT_PREFIX.counts.txt
$BEDTOOLS pairtopair -type both -a $OUTPUT_PREFIX.$REFNAME.segDup.PASS.bedpe -b $BEDPE -is -slop 50 > $OUTPUT_PREFIX.compared.$REFNAME.segDup.PASS.bedpe
echo "pairtopair for $REFNAME segDup PASS vs combined calls, 50 bp slop" >> $OUTPUT_PREFIX.counts.txt
cut -f 1-6 $OUTPUT_PREFIX.compared.$REFNAME.segDup.PASS.bedpe | sort -u | wc -l >> $OUTPUT_PREFIX.counts.txt

#More lenient sv comparison
/bin/bash $BEDPETOBED $OUTPUT_PREFIX.$REFNAME.nonRep.bedpe $OUTPUT_PREFIX.$REFNAME.nonRep.bed
/bin/bash $BEDPETOBED $OUTPUT_PREFIX.$REFNAME.nonRep.PASS.bedpe $OUTPUT_PREFIX.$REFNAME.nonRep.PASS.bed
/bin/bash $BEDPETOBED $OUTPUT_PREFIX.$REFNAME.str.bedpe $OUTPUT_PREFIX.$REFNAME.str.bed
/bin/bash $BEDPETOBED $OUTPUT_PREFIX.$REFNAME.str.PASS.bedpe $OUTPUT_PREFIX.$REFNAME.str.PASS.bed
/bin/bash $BEDPETOBED $OUTPUT_PREFIX.$REFNAME.segDup.bedpe $OUTPUT_PREFIX.$REFNAME.segDup.bed
/bin/bash $BEDPETOBED $OUTPUT_PREFIX.$REFNAME.segDup.PASS.bedpe $OUTPUT_PREFIX.$REFNAME.segDup.PASS.bed

/bin/bash $BEDPETOBED $BEDPE $BEDPE.bed
echo "total $REFNAME nonRep bed" >> $OUTPUT_PREFIX.counts.txt
cut -f 1-3 $OUTPUT_PREFIX.$REFNAME.nonRep.bed | sort -u | wc -l >> $OUTPUT_PREFIX.counts.txt
echo "total $REFNAME str bed" >> $OUTPUT_PREFIX.counts.txt
cut -f 1-3 $OUTPUT_PREFIX.$REFNAME.str.bed | sort -u | wc -l >> $OUTPUT_PREFIX.counts.txt
echo "total $REFNAME segDup bed" >> $OUTPUT_PREFIX.counts.txt
cut -f 1-3 $OUTPUT_PREFIX.$REFNAME.segDup.bed | sort -u | wc -l >> $OUTPUT_PREFIX.counts.txt
echo "total $REFNAME nonRep PASS bed" >> $OUTPUT_PREFIX.counts.txt
cut -f 1-3 $OUTPUT_PREFIX.$REFNAME.nonRep.PASS.bed | sort -u | wc -l >> $OUTPUT_PREFIX.counts.txt
echo "total $REFNAME str PASS bed" >> $OUTPUT_PREFIX.counts.txt
cut -f 1-3 $OUTPUT_PREFIX.$REFNAME.str.PASS.bed | sort -u | wc -l >> $OUTPUT_PREFIX.counts.txt
echo "total $REFNAME segDup PASS bed" >> $OUTPUT_PREFIX.counts.txt
cut -f 1-3 $OUTPUT_PREFIX.$REFNAME.segDup.PASS.bed | sort -u | wc -l >> $OUTPUT_PREFIX.counts.txt

$BEDTOOLS intersect -wa -u -f .5 -r -a $OUTPUT_PREFIX.$REFNAME.nonRep.bed -b $BEDPE.bed > $OUTPUT_PREFIX.intersect.50.$REFNAME.nonRep.bed
$BEDTOOLS intersect -wa -u -f .5 -r -a $OUTPUT_PREFIX.$REFNAME.str.bed -b $BEDPE.bed > $OUTPUT_PREFIX.intersect.50.$REFNAME.str.bed
$BEDTOOLS intersect -wa -u -f .5 -r -a $OUTPUT_PREFIX.$REFNAME.segDup.bed -b $BEDPE.bed > $OUTPUT_PREFIX.intersect.50.$REFNAME.segDup.bed
$BEDTOOLS intersect -wa -u -f .5 -r -a $OUTPUT_PREFIX.$REFNAME.nonRep.PASS.bed -b $BEDPE.bed > $OUTPUT_PREFIX.intersect.50.$REFNAME.nonRep.PASS.bed
$BEDTOOLS intersect -wa -u -f .5 -r -a $OUTPUT_PREFIX.$REFNAME.str.PASS.bed -b $BEDPE.bed > $OUTPUT_PREFIX.intersect.50.$REFNAME.str.PASS.bed
$BEDTOOLS intersect -wa -u -f .5 -r -a $OUTPUT_PREFIX.$REFNAME.segDup.PASS.bed -b $BEDPE.bed > $OUTPUT_PREFIX.intersect.50.$REFNAME.segDup.PASS.bed

$BEDTOOLS intersect -wa -u -f .1 -r -a $OUTPUT_PREFIX.$REFNAME.nonRep.bed -b $BEDPE.bed > $OUTPUT_PREFIX.intersect.10.$REFNAME.nonRep.bed
$BEDTOOLS intersect -wa -u -f .1 -r -a $OUTPUT_PREFIX.$REFNAME.str.bed -b $BEDPE.bed > $OUTPUT_PREFIX.intersect.10.$REFNAME.str.bed
$BEDTOOLS intersect -wa -u -f .1 -r -a $OUTPUT_PREFIX.$REFNAME.segDup.bed -b $BEDPE.bed > $OUTPUT_PREFIX.intersect.10.$REFNAME.segDup.bed
$BEDTOOLS intersect -wa -u -f .1 -r -a $OUTPUT_PREFIX.$REFNAME.nonRep.PASS.bed -b $BEDPE.bed > $OUTPUT_PREFIX.intersect.10.$REFNAME.nonRep.PASS.bed
$BEDTOOLS intersect -wa -u -f .1 -r -a $OUTPUT_PREFIX.$REFNAME.str.PASS.bed -b $BEDPE.bed > $OUTPUT_PREFIX.intersect.10.$REFNAME.str.PASS.bed
$BEDTOOLS intersect -wa -u -f .1 -r -a $OUTPUT_PREFIX.$REFNAME.segDup.PASS.bed -b $BEDPE.bed > $OUTPUT_PREFIX.intersect.10.$REFNAME.segDup.PASS.bed

echo "Reciprocal overlap 50% $REFNAME nonRep" >> $OUTPUT_PREFIX.counts.txt
cut -f 1-3 $OUTPUT_PREFIX.intersect.50.$REFNAME.nonRep.bed | sort -u | wc -l >> $OUTPUT_PREFIX.counts.txt
echo "Reciprocal overlap 50% $REFNAME str" >> $OUTPUT_PREFIX.counts.txt
cut -f 1-3 $OUTPUT_PREFIX.intersect.50.$REFNAME.str.bed | sort -u | wc -l >> $OUTPUT_PREFIX.counts.txt
echo "Reciprocal overlap 50% $REFNAME segDup" >> $OUTPUT_PREFIX.counts.txt
cut -f 1-3 $OUTPUT_PREFIX.intersect.50.$REFNAME.segDup.bed | sort -u | wc -l >> $OUTPUT_PREFIX.counts.txt
echo "Reciprocal overlap 50% $REFNAME nonRep PASS" >> $OUTPUT_PREFIX.counts.txt
cut -f 1-3 $OUTPUT_PREFIX.intersect.50.$REFNAME.nonRep.PASS.bed | sort -u | wc -l >> $OUTPUT_PREFIX.counts.txt
echo "Reciprocal overlap 50% $REFNAME str PASS" >> $OUTPUT_PREFIX.counts.txt
cut -f 1-3 $OUTPUT_PREFIX.intersect.50.$REFNAME.str.PASS.bed | sort -u | wc -l >> $OUTPUT_PREFIX.counts.txt
echo "Reciprocal overlap 50% $REFNAME segDup PASS" >> $OUTPUT_PREFIX.counts.txt
cut -f 1-3 $OUTPUT_PREFIX.intersect.50.$REFNAME.segDup.PASS.bed | sort -u | wc -l >> $OUTPUT_PREFIX.counts.txt

echo "Reciprocal overlap 10% $REFNAME nonRep" >> $OUTPUT_PREFIX.counts.txt
cut -f 1-3 $OUTPUT_PREFIX.intersect.10.$REFNAME.nonRep.bed | sort -u | wc -l >> $OUTPUT_PREFIX.counts.txt
echo "Reciprocal overlap 10% $REFNAME str" >> $OUTPUT_PREFIX.counts.txt
cut -f 1-3 $OUTPUT_PREFIX.intersect.10.$REFNAME.str.bed | sort -u | wc -l >> $OUTPUT_PREFIX.counts.txt
echo "Reciprocal overlap 10% $REFNAME segDup" >> $OUTPUT_PREFIX.counts.txt
cut -f 1-3 $OUTPUT_PREFIX.intersect.10.$REFNAME.segDup.bed | sort -u | wc -l >> $OUTPUT_PREFIX.counts.txt
echo "Reciprocal overlap 10% $REFNAME nonRep PASS" >> $OUTPUT_PREFIX.counts.txt
cut -f 1-3 $OUTPUT_PREFIX.intersect.10.$REFNAME.nonRep.PASS.bed | sort -u | wc -l >> $OUTPUT_PREFIX.counts.txt
echo "Reciprocal overlap 10% $REFNAME str PASS" >> $OUTPUT_PREFIX.counts.txt
cut -f 1-3 $OUTPUT_PREFIX.intersect.10.$REFNAME.str.PASS.bed | sort -u | wc -l >> $OUTPUT_PREFIX.counts.txt
echo "Reciprocal overlap 10% $REFNAME segDup PASS" >> $OUTPUT_PREFIX.counts.txt
cut -f 1-3 $OUTPUT_PREFIX.intersect.10.$REFNAME.segDup.PASS.bed | sort -u | wc -l >> $OUTPUT_PREFIX.counts.txt

echo "SVs in confident region intersect" >> $OUTPUT_PREFIX.counts.txt
$BEDTOOLS intersect -wa -u -f 1 -a $OUTPUT_PREFIX.ours.bed -b $REF_SV_BED | cut -f 1-3 | sort -u | wc -l >> $OUTPUT_PREFIX.counts.txt
echo "SVs in confident region both ends" >> $OUTPUT_PREFIX.counts.txt
$BEDTOOLS pairtobed -a $BEDPE -b $REF_SV_BED -type both | cut -f 1-6 | sort -u | wc -l >> $OUTPUT_PREFIX.counts.txt
echo "SVs in confident region coverage gt .9" >> $OUTPUT_PREFIX.counts.txt
$BEDTOOLS coverage -a $OUTPUT_PREFIX.ours.bed -b $REF_SV_BED | awk '$7 >= .9' | wc -l >> $OUTPUT_PREFIX.counts.txt

STR_BEDPE=$OUTPUT_PREFIX.ours.str.bedpe
SEGDUP_BEDPE=$OUTPUT_PREFIX.ours.segDup.bedpe
$BEDTOOLS pairtobed -a $BEDPE -b $STR_BED -type either > $STR_BEDPE
$BEDTOOLS pairtobed -a $BEDPE -b $SEGDUP_BED -type either > $SEGDUP_BEDPE
$PERL $COUNT_REPETITIVE $OUTPUT_PREFIX.$REFNAME.str.PASS.bedpe $OUTPUT_PREFIX.$REFNAME.str.bedpe $STR_BEDPE > $OUTPUT_PREFIX.str.VennInput.txt
$PERL $COUNT_REPETITIVE $OUTPUT_PREFIX.$REFNAME.segDup.PASS.bedpe $OUTPUT_PREFIX.$REFNAME.segDup.bedpe $SEGDUP_BEDPE > $OUTPUT_PREFIX.segDup.VennInput.txt

cat $OUTPUT_PREFIX.counts.txt | paste - - > $OUTPUT_PREFIX.counts.horizontal.txt
