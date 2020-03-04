## Assembly evaluation pipeline
This pipeline evaluates assemblies by aligning them to a reference, calling variants, quantifying those variants, and, if applicable, comparing them to gold-standard variant lists.  The end result is a set of plots.

### Inputs
* assembly_list - tab-delimited file with the fields <assembly name> <first fasta> <second fasta> <truth name>  If an assembly contains only one fasta, the second fasta can be set to an empty file.  If no truthset is available, <truth name> should be NONE.  Note that the <assembly name> should consist of the sample name followed by a period, followed by arbitrary other text.  Sample name should not contain a period.
* populations - tab-delimited file connecting sample to population.  Fields are: <sample> <population code> <superpopulation code> <sex>
* ref - reference fasta (not gzipped!)
* ref_index - fai index of reference fasta
* truth_list_snp_indel - tab-delimited file with the fields <truth name> <vcf location>
* truth_list_confident_regions - tab-delimited file with the fields <truth name> <bed location>
* truth_list_sv_vcf - tab-delimited file with the fields <truth name> <vcf location>
* gap_file - bed file denoting reference gaps
* str_track - bed file denoting simple repeats
* seg_dup_track - bed file denoting segmental duplications

### Alignment and variant calling methodology


### Outputs
* happy.sensitivity.png - Gives sensitivity to the snp/indel truthset as determined by hap.py
* het_fates.png - For each variant in the truthset, show whether it was called heterozygous, homozygous alt, or not called (homozygous ref) in the variant callset from the assembly.
* <assembly name>.no_gaps_coverage_summary.png - Gives the percentage of chr1-22 non-gap reference bases at each coverage 0X-9X, as well as 10X+
* <assembly name>.no_gaps_coverage.png - Same thing, broken down by chr and including chrX and chrY
* <assembly name>.ours.nonRep.lengths.txt.sv_lengths.png - Histogram of INS and DEL lengths.  Between 0 and 1000 bp the bin size is 10bp.  There are also bins for 1001-10000bp, 10001-100000bp, 100001-1000000bp, and 1Mb+.
* nonrepSvCounts.png - Counts for all SVs that weren't classified as STR or SegDup, broken down by SV type
* segDupSvCounts.png - Counts for all SVs that weren't classified as STR but did match the SegDup track with either breakpoint.
* strSvCounts.png - Counts for all SVs that matched the STR track with either breakpoint.
* sv.sensitivity.png - Percentage of the SV truthset found in the assembly using a variety of comparison techniques.  SVs are stratified by whether either breakpoint matched the STR track, SegDup track, or neither (Non-repetitive)
* sv.sensitivityCounts.png - Gives the counts of SVs from the truthset detected by each comparison method.
* truthHetIndelFates.png - For each heterozygous indel in the truthset, show whether it was called heterozygous, homozygous alt, or not called (homozygous ref) using 1bp, 10bp, or 50bp slop.
* truthHetIndelFatesErrorsOnly.png - Same as above but only shows the errors to give more detailed counts.
* truthset.repetitive.count.png - Counts the number of unique elements from the STR and SegDup tracks covered by the SV truthset, and determines how many of those were also found in the SV callset.
* truthset.repetitive.sensitivity.png - Same as above but displayed as sensitivity.
