#!/bin/bash

REF_INDEL_VCF=$1
INDEL_VCF=$2
OUTPUT_PREFIX=$3
SAMPLE=$4
REFNAME=$5
SEGDUP_BED=$6
STR_BED=$7

BEDTOOLS=/opt/hall-lab/bedtools
GREP=/bin/grep
PYTHON=/opt/hall-lab/python-2.7.15/bin/python
VCFTOBEDPE=/opt/hall-lab/scripts/vcfToBedpe.py
#SEGDUP_BED=/gscmnt/gc2802/halllab/aregier/jira/RI-477/v2/align_contigs_to_ref/HG002.GRCh38_noalts.segDup.bed
#STR_BED=/gscmnt/gc2802/halllab/aregier/jira/RI-477/v2/align_contigs_to_ref/HG002.GRCh38_noalts.str.merged.bed

#Preprocess ref vcf
$PYTHON $VCFTOBEDPE -i $REF_INDEL_VCF -o $OUTPUT_PREFIX.$REFNAME.indels.bedpe -m 1
$BEDTOOLS pairtobed -a $OUTPUT_PREFIX.$REFNAME.indels.bedpe -b $STR_BED -type either > $OUTPUT_PREFIX.$REFNAME.indels.str.bedpe
$BEDTOOLS pairtobed -a $OUTPUT_PREFIX.$REFNAME.indels.bedpe -b $SEGDUP_BED -type either > $OUTPUT_PREFIX.$REFNAME.allSegDup.bedpe
$BEDTOOLS pairtobed -a <($BEDTOOLS pairtobed -a $OUTPUT_PREFIX.$REFNAME.indels.bedpe -b $STR_BED -type neither) -b $SEGDUP_BED -type either > $OUTPUT_PREFIX.$REFNAME.indels.segdup.bedpe
$BEDTOOLS pairtobed -a <($BEDTOOLS pairtobed -a $OUTPUT_PREFIX.$REFNAME.indels.bedpe -b $STR_BED -type neither) -b $SEGDUP_BED -type neither > $OUTPUT_PREFIX.$REFNAME.indels.nonRep.bedpe

#Preprocess query vcf
$PYTHON $VCFTOBEDPE -i $INDEL_VCF -o $OUTPUT_PREFIX.indels.bedpe -m 1

echo "total $REFNAME indels" > $OUTPUT_PREFIX.indel.counts.txt
$GREP -v "^#" $OUTPUT_PREFIX.$REFNAME.indels.bedpe | sort -u | wc -l >> $OUTPUT_PREFIX.indel.counts.txt
echo "total $REFNAME het indels" > $OUTPUT_PREFIX.indel.counts.txt
$GREP -v "^#" $OUTPUT_PREFIX.$REFNAME.indels.bedpe | $GREP HET | sort -u | wc -l >> $OUTPUT_PREFIX.indel.counts.txt

$BEDTOOLS pairtopair -type both -a $OUTPUT_PREFIX.$REFNAME.indels.bedpe -b $OUTPUT_PREFIX.indels.bedpe -is -slop 50 > $OUTPUT_PREFIX.compared.$REFNAME.indels.50.bedpe
echo "pairtopair for $REFNAME indels vs our calls, 50 bp slop" >> $OUTPUT_PREFIX.indel.counts.txt
cut -f 1-6 $OUTPUT_PREFIX.compared.$REFNAME.indels.50.bedpe | sort -u | wc -l >> $OUTPUT_PREFIX.indel.counts.txt
echo "pairtopair for $REFNAME het indels vs our calls, 50 bp slop" >> $OUTPUT_PREFIX.indel.counts.txt
cut -f 1-7 $OUTPUT_PREFIX.compared.$REFNAME.indels.50.bedpe | $GREP HET | cut -f 1-6 | sort -u | wc -l >> $OUTPUT_PREFIX.indel.counts.txt
echo "pairtopair for $REFNAME het to het, 50 bp slop" >> $OUTPUT_PREFIX.indel.counts.txt
$GREP "HET.*HET" $OUTPUT_PREFIX.compared.$REFNAME.indels.50.bedpe | cut -f 1-6 | sort -u | wc -l >> $OUTPUT_PREFIX.indel.counts.txt
echo "pairtopair for $REFNAME het to homalt, 50 bp slop" >> $OUTPUT_PREFIX.indel.counts.txt
$GREP "HET.*HOMALT" $OUTPUT_PREFIX.compared.$REFNAME.indels.50.bedpe | cut -f 1-6 | sort -u | wc -l >> $OUTPUT_PREFIX.indel.counts.txt

$BEDTOOLS pairtopair -type both -a $OUTPUT_PREFIX.$REFNAME.indels.bedpe -b $OUTPUT_PREFIX.indels.bedpe -is -slop 10 > $OUTPUT_PREFIX.compared.$REFNAME.indels.10.bedpe
echo "pairtopair for $REFNAME indels vs our calls, 10 bp slop" >> $OUTPUT_PREFIX.indel.counts.txt
cut -f 1-6 $OUTPUT_PREFIX.compared.$REFNAME.indels.10.bedpe | sort -u | wc -l >> $OUTPUT_PREFIX.indel.counts.txt
echo "pairtopair for $REFNAME het indels vs our calls, 10 bp slop" >> $OUTPUT_PREFIX.indel.counts.txt
cut -f 1-7 $OUTPUT_PREFIX.compared.$REFNAME.indels.10.bedpe | grep HET | cut -f 1-6 | sort -u | wc -l >> $OUTPUT_PREFIX.indel.counts.txt
echo "pairtopair for $REFNAME het to het, 10 bp slop" >> $OUTPUT_PREFIX.indel.counts.txt
$GREP "HET.*HET" $OUTPUT_PREFIX.compared.$REFNAME.indels.10.bedpe | cut -f 1-6 | sort -u | wc -l >> $OUTPUT_PREFIX.indel.counts.txt
echo "pairtopair for $REFNAME het to homalt, 10 bp slop" >> $OUTPUT_PREFIX.indel.counts.txt
$GREP "HET.*HOMALT" $OUTPUT_PREFIX.compared.$REFNAME.indels.10.bedpe | cut -f 1-6 | sort -u | wc -l >> $OUTPUT_PREFIX.indel.counts.txt

$BEDTOOLS pairtopair -type both -a $OUTPUT_PREFIX.$REFNAME.indels.bedpe -b $OUTPUT_PREFIX.indels.bedpe -is -slop 1 > $OUTPUT_PREFIX.compared.$REFNAME.indels.1.bedpe
echo "pairtopair for $REFNAME indels vs our calls, 1 bp slop" >> $OUTPUT_PREFIX.indel.counts.txt
cut -f 1-6 $OUTPUT_PREFIX.compared.$REFNAME.indels.1.bedpe | sort -u | wc -l >> $OUTPUT_PREFIX.indel.counts.txt
echo "pairtopair for $REFNAME het indels vs our calls, 1 bp slop" >> $OUTPUT_PREFIX.indel.counts.txt
cut -f 1-7 $OUTPUT_PREFIX.compared.$REFNAME.indels.1.bedpe | $GREP HET | cut -f 1-6 | sort -u | wc -l >> $OUTPUT_PREFIX.indel.counts.txt
echo "pairtopair for $REFNAME het to het, 1 bp slop" >> $OUTPUT_PREFIX.indel.counts.txt
$GREP "HET.*HET" $OUTPUT_PREFIX.compared.$REFNAME.indels.1.bedpe | cut -f 1-6 | sort -u | wc -l >> $OUTPUT_PREFIX.indel.counts.txt
echo "pairtopair for $REFNAME het to homalt, 1 bp slop" >> $OUTPUT_PREFIX.indel.counts.txt
$GREP "HET.*HOMALT" $OUTPUT_PREFIX.compared.$REFNAME.indels.1.bedpe | cut -f 1-6 | sort -u | wc -l >> $OUTPUT_PREFIX.indel.counts.txt
cat $OUTPUT_PREFIX.indel.counts.txt | paste - - > $OUTPUT_PREFIX.indel.counts.horizontal.txt
