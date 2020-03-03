#!/bin/bash
set -exo pipefail
CONTIG1=$1
CONTIG2=$2
REF=$3
OUTPUT_PREFIX=$4
SAMPLE=$5
GAP_FILE=$6
STR_TRACK=$7
SEG_DUP_TRACK=$8

MINIMAP2=/opt/hall-lab/minimap2/minimap2
K8=/opt/hall-lab/minimap2/k8
PAFTOOLS=/opt/hall-lab/minimap2/misc/paftools.js
SAMTOOLS=/opt/hall-lab/samtools-1.9/bin/samtools
PYTHON=/opt/hall-lab/python-2.7.15/bin/python
VAR_TO_VCF=/opt/hall-lab/scripts/varToVcf.py
VAR_TO_BEDPE=/opt/hall-lab/scripts/varToBedpe.py
SPLIT_TO_BEDPE=/opt/hall-lab/scripts/splitReadSamToBedpe
BEDPE_TO_BKPTS=/opt/hall-lab/scripts/splitterToBreakpoint
BEDPETOBED=/opt/hall-lab/scripts/bedpetobed.sh
REARRANGE_BREAKPOINTS=/opt/hall-lab/scripts/rearrange_breakpoints.pl
ADD_ALIGNMENT_GAP_INFO=/opt/hall-lab/scripts/add_alignment_gap_info.pl
SUMMARIZE_GENOME_COVERAGE=/opt/hall-lab/scripts/summarize_genome_coverage.py
SVTOOLS=/opt/hall-lab/python-2.7.15/bin/svtools
PERL=/usr/bin/perl
BGZIP=/opt/hall-lab/htslib-1.9/bin/bgzip
TABIX=/opt/hall-lab/htslib-1.9/bin/tabix
GENOTYPE_VCF=/opt/hall-lab/scripts/vcfToGenotyped.pl
BEDTOOLS=/opt/hall-lab/bedtools
GREP=/bin/grep
#GAP_FILE=/gscmnt/gc2802/halllab/ccdg_resources/genomes/human/GRCh38DH/annotations/ira_b38_061117/download_061117/gap.053117.b38.sorted.bed
#STR_TRACK=/gscmnt/gc2802/halllab/ccdg_resources/genomes/human/GRCh38DH/annotations/ira_b38_061117/download_061117/simpleRepeat.b38.sorted.bed
#SEG_DUP_TRACK=/gscmnt/gc2802/halllab/ccdg_resources/genomes/human/GRCh38DH/annotations/ira_b38_061117/download_061117/genomicSuperDups.b38.sorted.bed

CONTIGS=$OUTPUT_PREFIX.combined_contigs.fa
if [ -s $CONTIG2 ]; then
    cat <(cat $CONTIG1 | sed 's/\(>.*\)/\1_H1/') <(cat $CONTIG2 | sed 's/\(>.*\)/\1_H2/') > $CONTIGS
else
    cp $CONTIG1 $CONTIGS
fi
$MINIMAP2 -x asm5 --cs $REF $CONTIGS > $OUTPUT_PREFIX.paf
$MINIMAP2 -ax asm5 -L --cs $REF $CONTIGS | $SAMTOOLS sort -T $OUTPUT_PREFIX.tmp -O bam - > $OUTPUT_PREFIX.bam
sort -k6,6 -k8,8n $OUTPUT_PREFIX.paf | $K8 $PAFTOOLS call -f $REF -s $SAMPLE -l 1 -L 1 -q 0 - | $BGZIP -c > $OUTPUT_PREFIX.loose.var.vcf.gz
$TABIX -fp vcf $OUTPUT_PREFIX.loose.var.vcf.gz
sort -k6,6 -k8,8n $OUTPUT_PREFIX.paf | $K8 $PAFTOOLS call -l 1 -L 1 -q 0 - | grep "^V" | sort -V | $BGZIP -c > $OUTPUT_PREFIX.loose.var.txt.gz
$PYTHON $VAR_TO_VCF -i <(zcat $OUTPUT_PREFIX.loose.var.txt.gz) -r $REF -s $SAMPLE -o $OUTPUT_PREFIX.loose.vcf
$SVTOOLS vcfsort $OUTPUT_PREFIX.loose.vcf | $PERL $GENOTYPE_VCF | $BGZIP -c > $OUTPUT_PREFIX.loose.genotyped.vcf.gz
$TABIX -f -p vcf $OUTPUT_PREFIX.loose.genotyped.vcf.gz

#Get putative SV breakpoints
$SAMTOOLS sort -n -T $OUTPUT_PREFIX.tmp -O bam $OUTPUT_PREFIX.bam > $OUTPUT_PREFIX.namesorted.bam
$SAMTOOLS view -h -F 4 $OUTPUT_PREFIX.namesorted.bam | $PYTHON $SPLIT_TO_BEDPE -i stdin > $OUTPUT_PREFIX.bedpe

$PYTHON $BEDPE_TO_BKPTS -i $OUTPUT_PREFIX.bedpe -f $OUTPUT_PREFIX -q $CONTIGS -e $REF > $OUTPUT_PREFIX.breakpoints.bedpe
$SVTOOLS bedpesort $OUTPUT_PREFIX.breakpoints.bedpe | $PERL $REARRANGE_BREAKPOINTS > $OUTPUT_PREFIX.breakpoints.sorted.bedpe
cat <($GREP "^#" $OUTPUT_PREFIX.breakpoints.sorted.bedpe) <(paste <($GREP -v "^#" $OUTPUT_PREFIX.breakpoints.sorted.bedpe | cut -f 1-6) <(paste -d : <($GREP -v "^#" $OUTPUT_PREFIX.breakpoints.sorted.bedpe | cut -f 7) <($GREP -v "^#" $OUTPUT_PREFIX.breakpoints.sorted.bedpe | cut -f 19 | sed 's/.*SVLEN=/SVLEN=/' | sed 's/;.*//')) <($GREP -v "^#" $OUTPUT_PREFIX.breakpoints.sorted.bedpe | cut -f 8-)) | $PERL $ADD_ALIGNMENT_GAP_INFO> $OUTPUT_PREFIX.breakpoints.sorted.fixed.bedpe
mv $OUTPUT_PREFIX.breakpoints.sorted.fixed.bedpe $OUTPUT_PREFIX.breakpoints.sorted.bedpe
$SVTOOLS bedpetobed12 -i $OUTPUT_PREFIX.breakpoints.sorted.bedpe -o $OUTPUT_PREFIX.breakpoints.bed -n $SAMPLE

#Alignment stats
$BEDTOOLS genomecov -ibam $OUTPUT_PREFIX.bam > $OUTPUT_PREFIX.genomecov.hist
$BEDTOOLS genomecov -ibam $OUTPUT_PREFIX.bam -bga > $OUTPUT_PREFIX.genomecov.bga
$BEDTOOLS intersect -v -a $OUTPUT_PREFIX.genomecov.bga -b <(grep -v "^track" $GAP_FILE | sed 's/^/chr/') | $PYTHON $SUMMARIZE_GENOME_COVERAGE > $OUTPUT_PREFIX.genomecov.noGaps.hist

#Counts
echo "Cigar substitutions" > $OUTPUT_PREFIX.counts.txt
zcat $OUTPUT_PREFIX.loose.genotyped.vcf.gz | $GREP -v "^#" | cut -f 1-5 | sort -u | cut -f 4-5 | $GREP '^.	.$' | wc -l >> $OUTPUT_PREFIX.counts.txt
echo "Cigar all deletions" >> $OUTPUT_PREFIX.counts.txt
zcat $OUTPUT_PREFIX.loose.genotyped.vcf.gz | $GREP -v "^#" | cut -f 1-5 | sort -u | cut -f 4-5 | $GREP -v '^.	.$' | $GREP '	.$' | wc -l >> $OUTPUT_PREFIX.counts.txt
echo "Cigar all insertions" >> $OUTPUT_PREFIX.counts.txt
zcat $OUTPUT_PREFIX.loose.genotyped.vcf.gz | $GREP -v "^#" | cut -f 1-5 | sort -u | cut -f 4-5 | $GREP -v '^.	.$' | $GREP '^.	' | wc -l >> $OUTPUT_PREFIX.counts.txt
echo "Cigar 1bp deletions" >> $OUTPUT_PREFIX.counts.txt
zcat $OUTPUT_PREFIX.loose.genotyped.vcf.gz | $GREP -v "^#" | cut -f 1-5 | sort -u | cut -f 4-5 | $GREP '^..	.$' | wc -l >> $OUTPUT_PREFIX.counts.txt
echo "Cigar 1bp insertions" >> $OUTPUT_PREFIX.counts.txt
zcat $OUTPUT_PREFIX.loose.genotyped.vcf.gz | $GREP -v "^#" | cut -f 1-5 | sort -u | cut -f 4-5 | $GREP '^.	..$' | wc -l >> $OUTPUT_PREFIX.counts.txt
echo "Cigar 2bp deletions" >> $OUTPUT_PREFIX.counts.txt
zcat $OUTPUT_PREFIX.loose.genotyped.vcf.gz | $GREP -v "^#" | cut -f 1-5 | sort -u | cut -f 4-5 | $GREP '^...	.$' | wc -l >> $OUTPUT_PREFIX.counts.txt
echo "Cigar 2bp insertions" >> $OUTPUT_PREFIX.counts.txt
zcat $OUTPUT_PREFIX.loose.genotyped.vcf.gz | $GREP -v "^#" | cut -f 1-5 | sort -u | cut -f 4-5 | $GREP '^.	...$' | wc -l >> $OUTPUT_PREFIX.counts.txt
echo "Split breakpoint dels" >> $OUTPUT_PREFIX.counts.txt
cat $OUTPUT_PREFIX.breakpoints.sorted.bedpe | $GREP "DEL" | cut -f 1-6 | sort -u | wc -l >> $OUTPUT_PREFIX.counts.txt
echo "Split breakpoint dups" >> $OUTPUT_PREFIX.counts.txt
cat $OUTPUT_PREFIX.breakpoints.sorted.bedpe | $GREP "DUP" | cut -f 1-6 | sort -u | wc -l >> $OUTPUT_PREFIX.counts.txt
echo "Split breakpoint invs" >> $OUTPUT_PREFIX.counts.txt
cat $OUTPUT_PREFIX.breakpoints.sorted.bedpe | $GREP "INV" | cut -f 1-6 | sort -u | wc -l >> $OUTPUT_PREFIX.counts.txt
#Combine bedpes
$PYTHON $VAR_TO_BEDPE -i <(zcat $OUTPUT_PREFIX.loose.var.txt.gz) -m 40 -s $SAMPLE -o $OUTPUT_PREFIX.loose.indels40plus.bedpe
$SVTOOLS bedpesort $OUTPUT_PREFIX.loose.indels40plus.bedpe | uniq > $OUTPUT_PREFIX.loose.indels40plus.sorted.bedpe
$SVTOOLS bedpetobed12 -i $OUTPUT_PREFIX.loose.indels40plus.sorted.bedpe -o $OUTPUT_PREFIX.loose.indels40plus.sorted.forIGV.bed
$SVTOOLS bedpesort <(cat <(cat $OUTPUT_PREFIX.breakpoints.sorted.bedpe | sort -u | cut -f 1-22) <($GREP -v "^#" $OUTPUT_PREFIX.loose.indels40plus.sorted.bedpe | sort -u)) $OUTPUT_PREFIX.ours.bedpe

$BEDTOOLS merge -i $STR_TRACK | sed 's/^/chr/'> $OUTPUT_PREFIX.str.merged.bed
$BEDTOOLS merge -i $SEG_DUP_TRACK | $GREP -v "^t" | sed 's/^/chr/' | cut -f 1-3 > $OUTPUT_PREFIX.segDup.bed

#Split out STR and segDup
$BEDTOOLS pairtobed -a $OUTPUT_PREFIX.ours.bedpe -b $OUTPUT_PREFIX.str.merged.bed -type either > $OUTPUT_PREFIX.ours.str.bedpe
$BEDTOOLS pairtobed -a $OUTPUT_PREFIX.ours.bedpe -b $OUTPUT_PREFIX.segDup.bed -type either > $OUTPUT_PREFIX.ours.allSegDup.bedpe
$BEDTOOLS pairtobed -a <($BEDTOOLS pairtobed -a $OUTPUT_PREFIX.ours.bedpe -b $OUTPUT_PREFIX.str.merged.bed -type neither) -b $OUTPUT_PREFIX.segDup.bed -type either > $OUTPUT_PREFIX.ours.segDup.bedpe
$BEDTOOLS pairtobed -a <($BEDTOOLS pairtobed -a $OUTPUT_PREFIX.ours.bedpe -b $OUTPUT_PREFIX.str.merged.bed -type neither) -b $OUTPUT_PREFIX.segDup.bed -type neither > $OUTPUT_PREFIX.ours.nonRep.bedpe
echo "total combined calls" >> $OUTPUT_PREFIX.counts.txt
$GREP -v "^#" $OUTPUT_PREFIX.ours.bedpe | cut -f 1-6 | sort -u | wc -l >> $OUTPUT_PREFIX.counts.txt
echo "total nonRep combined calls" >> $OUTPUT_PREFIX.counts.txt
$GREP -v "^#" $OUTPUT_PREFIX.ours.nonRep.bedpe | cut -f 1-6 | sort -u | wc -l >> $OUTPUT_PREFIX.counts.txt
echo "total str combined calls" >> $OUTPUT_PREFIX.counts.txt
$GREP -v "^#" $OUTPUT_PREFIX.ours.str.bedpe | cut -f 1-6 | sort -u | wc -l >> $OUTPUT_PREFIX.counts.txt
echo "total segDup combined calls" >> $OUTPUT_PREFIX.counts.txt
$GREP -v "^#" $OUTPUT_PREFIX.ours.segDup.bedpe | cut -f 1-6 | sort -u | wc -l >> $OUTPUT_PREFIX.counts.txt

cut -f 1-6,11,19 $OUTPUT_PREFIX.ours.nonRep.bedpe | sort -u | sed 's/SVTYPE=INS//' | sed 's/SVTYPE=DEL//' | sed 's/SVTYPE=TRANS//' | sed 's/SVTYPE=GAP//' | sed 's/DEL=[0-9]*;//' | sed 's/REF_OLAP=[0-9]*;//' | sed 's/OVERLAP=-[0-9]*//' | sed 's/OVERLAP=[0-9]*//' | sed 's/INS=[0-9]*//' | sed 's/GAP=[0-9]*//' | sed 's/SVLEN=//' | sed 's/SVTYPE=DUP//' | sed 's/SVTYPE=G_OVERLAP//' | sed 's/SVTYPE=INV//' | sed 's/REF_OLAP=-[0-9]*//' | sed 's/;//g' | cut -f 7-8 | grep -v TRANS | grep -v DUP | grep -v INV | grep -v G_OVERLAP > $OUTPUT_PREFIX.ours.nonRep.lengths.txt
cut -f 1-6,11,19 $OUTPUT_PREFIX.ours.str.bedpe | sort -u | sed 's/SVTYPE=INS//' | sed 's/SVTYPE=DEL//' | sed 's/SVTYPE=TRANS//' | sed 's/SVTYPE=GAP//' | sed 's/DEL=[0-9]*;//' | sed 's/REF_OLAP=[0-9]*;//' | sed 's/OVERLAP=-[0-9]*//' | sed 's/OVERLAP=[0-9]*//' | sed 's/INS=[0-9]*//' | sed 's/GAP=[0-9]*//' | sed 's/SVLEN=//' | sed 's/SVTYPE=DUP//' | sed 's/SVTYPE=G_OVERLAP//' | sed 's/SVTYPE=INV//' | sed 's/REF_OLAP=-[0-9]*//' | sed 's/;//g' | cut -f 7-8 | grep -v TRANS | grep -v DUP | grep -v INV | grep -v G_OVERLAP > $OUTPUT_PREFIX.ours.str.lengths.txt
cut -f 1-6,11,19 $OUTPUT_PREFIX.ours.segDup.bedpe | sort -u | sed 's/SVTYPE=INS//' | sed 's/SVTYPE=DEL//' | sed 's/SVTYPE=TRANS//' | sed 's/SVTYPE=GAP//' | sed 's/DEL=[0-9]*;//' | sed 's/REF_OLAP=[0-9]*;//' | sed 's/OVERLAP=-[0-9]*//' | sed 's/OVERLAP=[0-9]*//' | sed 's/INS=[0-9]*//' | sed 's/GAP=[0-9]*//' | sed 's/SVLEN=//' | sed 's/SVTYPE=DUP//' | sed 's/SVTYPE=G_OVERLAP//' | sed 's/SVTYPE=INV//' | sed 's/REF_OLAP=-[0-9]*//' | sed 's/;//g' | cut -f 7-8 | grep -v TRANS | grep -v DUP | grep -v INV | grep -v G_OVERLAP > $OUTPUT_PREFIX.ours.segDup.lengths.txt

/bin/bash $BEDPETOBED $OUTPUT_PREFIX.ours.nonRep.bedpe $OUTPUT_PREFIX.ours.del.bed
$BEDTOOLS sort -i $OUTPUT_PREFIX.ours.del.bed | uniq > $OUTPUT_PREFIX.ours.tmp
mv $OUTPUT_PREFIX.ours.tmp $OUTPUT_PREFIX.ours.del.bed
echo "total our calls combined deletions" >> $OUTPUT_PREFIX.counts.txt
cut -f 1-3 $OUTPUT_PREFIX.ours.del.bed | sort -u | wc -l >> $OUTPUT_PREFIX.counts.txt

$BEDTOOLS merge -i $OUTPUT_PREFIX.ours.del.bed > $OUTPUT_PREFIX.ours.del.merged.bed
echo "total our calls combined deletions merged" >> $OUTPUT_PREFIX.counts.txt
cut -f 1-3 $OUTPUT_PREFIX.ours.del.merged.bed | sort -u | wc -l  >> $OUTPUT_PREFIX.counts.txt

$GREP GAP_DIFF $OUTPUT_PREFIX.breakpoints.sorted.bedpe | $GREP INS | sed 's/\t[a-zA-Z0-9_\-]*:/\t/' | cut -f 1-7,11,19 | sed 's/INS=[0-9]*//' | sort -u | sed 's/.*\tSVLEN=//' | sed 's/:GAP_DIFF=/\t/' | sed 's/:R_OLAP=/\t/' | sed 's/:QUERY_OLAP=/\t/' | sed 's/INS\t;REF_OLAP=//' | sed 's/;.*//' > $OUTPUT_PREFIX.breakpoints.INS.gapRatio.txt
$GREP GAP_DIFF $OUTPUT_PREFIX.breakpoints.sorted.bedpe | $GREP DEL | sed 's/\t[a-zA-Z0-9_\-]*:/\t/' | cut -f 1-7,11,19 | sed 's/DEL=[0-9]*//' | sort -u | sed 's/.*\tSVLEN=//' | sed 's/:GAP_DIFF=/\t/' | sed 's/:R_OLAP=/\t/' | sed 's/:QUERY_OLAP=/\t/' | sed 's/DEL\t;REF_OLAP=//' | sed 's/;.*//' > $OUTPUT_PREFIX.breakpoints.DEL.gapRatio.txt

cut -f 1-6,11 $OUTPUT_PREFIX.ours.nonRep.bedpe | sort -u | cut -f 7 | sort | uniq -c | sed 's/^ *//' | sed 's/ +/\t/' > $OUTPUT_PREFIX.ours.nonRep.typeCount.txt
cut -f 1-6,11 $OUTPUT_PREFIX.ours.str.bedpe | sort -u | cut -f 7 | sort | uniq -c | sed 's/^ *//' | sed 's/ +/\t/' > $OUTPUT_PREFIX.ours.str.typeCount.txt
cut -f 1-6,11 $OUTPUT_PREFIX.ours.segDup.bedpe | sort -u | cut -f 7 | sort | uniq -c | sed 's/^ *//' | sed 's/ +/\t/' > $OUTPUT_PREFIX.ours.segDup.typeCount.txt

cat $OUTPUT_PREFIX.counts.txt | paste - - > $OUTPUT_PREFIX.counts.horizontal.txt
