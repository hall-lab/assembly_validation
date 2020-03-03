#!/bin/bash

VCF=$1
OUTPUT_PREFIX=$2

echo "SNP	HET	HET" > $OUTPUT_PREFIX.het.counts.txt
zcat $VCF | grep -v "^#" | grep SNP:het.*SNP:het | wc -l >> $OUTPUT_PREFIX.het.counts.txt
echo "SNP	HET	HOMREF" >> $OUTPUT_PREFIX.het.counts.txt
zcat $VCF | grep -v "^#" | grep SNP:het.*NO | wc -l >> $OUTPUT_PREFIX.het.counts.txt
echo "SNP	HET	HOMALT" >> $OUTPUT_PREFIX.het.counts.txt
zcat $VCF | grep -v "^#" | grep SNP:het.*SNP:homalt | wc -l >> $OUTPUT_PREFIX.het.counts.txt
echo "INDEL	HET	HET" >> $OUTPUT_PREFIX.het.counts.txt
zcat $VCF | grep -v "^#" | grep INDEL:het.*INDEL:het | wc -l >> $OUTPUT_PREFIX.het.counts.txt
echo "INDEL	HET	HOMREF" >> $OUTPUT_PREFIX.het.counts.txt
zcat $VCF | grep -v "^#" | grep INDEL:het.*NO | wc -l >> $OUTPUT_PREFIX.het.counts.txt
echo "INDEL	HET	HOMALT" >> $OUTPUT_PREFIX.het.counts.txt
zcat $VCF | grep -v "^#" | grep INDEL:het.*INDEL:homalt | wc -l >> $OUTPUT_PREFIX.het.counts.txt
echo "SNP	HOMALT	HET" >> $OUTPUT_PREFIX.het.counts.txt
zcat $VCF | grep -v "^#" | grep SNP:hom.*SNP:het | wc -l >> $OUTPUT_PREFIX.het.counts.txt
echo "SNP	HOMALT	HOMREF" >> $OUTPUT_PREFIX.het.counts.txt
zcat $VCF | grep -v "^#" | grep SNP:hom.*NO | wc -l >> $OUTPUT_PREFIX.het.counts.txt
echo "SNP	HOMALT	HOMALT" >> $OUTPUT_PREFIX.het.counts.txt
zcat $VCF | grep -v "^#" | grep SNP:hom.*SNP:homalt | wc -l >> $OUTPUT_PREFIX.het.counts.txt
echo "INDEL	HOMALT	HET" >> $OUTPUT_PREFIX.het.counts.txt
zcat $VCF | grep -v "^#" | grep INDEL:hom.*INDEL:het | wc -l >> $OUTPUT_PREFIX.het.counts.txt
echo "INDEL	HOMALT	HOMREF" >> $OUTPUT_PREFIX.het.counts.txt
zcat $VCF | grep -v "^#" | grep INDEL:hom.*NO | wc -l >> $OUTPUT_PREFIX.het.counts.txt
echo "INDEL	HOMALT	HOMALT" >> $OUTPUT_PREFIX.het.counts.txt
zcat $VCF | grep -v "^#" | grep INDEL:hom.*INDEL:homalt | wc -l >> $OUTPUT_PREFIX.het.counts.txt

cat $OUTPUT_PREFIX.het.counts.txt | paste - - > $OUTPUT_PREFIX.het.counts.horizontal.txt
