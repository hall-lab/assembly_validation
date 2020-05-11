library(data.table)
library(reshape2)
library(plyr)
library(ggplot2)
library(dplyr)
library(stringr)

args <- commandArgs(trailingOnly = TRUE)
populations <- read.table(file=args[1], header=FALSE, sep="\t", col.names=c("Sample", "Population", "Superpopulation", "Sex"))
small_variants_files <- read.table(file=args[2], header=FALSE, col.names=c("file"))
indels_files <- read.table(file=args[3], header=FALSE, col.names=c("file"))
sv_files <- read.table(file=args[4], header=FALSE, col.names=c("file"))
genome_cov_files <- read.table(file=args[5], header=FALSE, col.names=c("file"))
nonRep_lengths_files <- read.table(file=args[6], header=FALSE, col.names=c("file"))
str_lengths_files <- read.table(file=args[7], header=FALSE, col.names=c("file"))
segDup_lengths_files <- read.table(file=args[8], header=FALSE, col.names=c("file"))
counts_files <- read.table(file=args[9], header=FALSE, col.names=c("file"))
nonRep_types_files <- read.table(file=args[10], header=FALSE, col.names=c("file"))
str_types_files <- read.table(file=args[11], header=FALSE, col.names=c("file"))
segDup_types_files <- read.table(file=args[12], header=FALSE, col.names=c("file"))
het_fates_files <- read.table(file=args[13], header=FALSE, col.names=c("file"))
str_venn_files <- read.table(file=args[14], header=FALSE, col.names=c("file"))
segDup_venn_files <- read.table(file=args[15], header=FALSE, col.names=c("file"))
cigar_indel_lengths_file <- read.table(file=args[16], header=FALSE, col.names=c("file"))
assembly_file <- read.table(file=args[17], header=FALSE, col.names=c("Name", "Fasta1", "Fasta2", "Truth", "ID", "DataType"), sep="\t")
assemblies <- assembly_file[,c("Name", "ID", "DataType")]
assemblies$Squashed <- rep(FALSE, length(assemblies$Name))
assemblies$Squashed[grepl("empty.fa", assembly_file$Fasta2)] <- TRUE

cigar_indel_lengths <- data.frame()
for (in_file in cigar_indel_lengths_file$file) {
    assembly <- gsub(".indel_lengths.horizontal.txt", "", basename(in_file))
    df <- read.table(file=in_file, header=TRUE, sep="\t")
    df$Assembly <- assembly
    cigar_indel_lengths <- rbind(cigar_indel_lengths, df)
}
cigar_indel_lengths$Metric = factor(cigar_indel_lengths$Metric, levels=c("1 bp", "2 bp", "3-10 bp", "11-100 bp", "101-1000 bp", "1-10 kb", "10 kb+"))
head(assemblies)
cigar_indel_lengths <- merge(cigar_indel_lengths, assemblies, by.x="Assembly", by.y="Name")
ggplot(cigar_indel_lengths, aes(x=ID, y=Value, fill=Squashed)) +
    geom_bar(stat="identity", position="dodge") +
    facet_grid(Metric~Type, scales="free") +
    theme_bw() +
    theme(axis.text.x = element_text(angle=45, hjust=1, size=7)) +
    ggtitle("Counts of indel variants detected in cigar string by size")
ggsave("plots/indel_cigar_counts.png")

counts <- data.frame()
for (in_file in counts_files$file) {
    assembly <- gsub(".counts.horizontal.txt", "", basename(in_file))
    df <- read.table(file=in_file, header=FALSE, sep="\t", col.names=c("Metric", "Value"))
    df$Assembly <- assembly
    counts <- rbind(counts, df)
}

counts <- merge(counts, assemblies, by.x="Assembly", by.y="Name")
ggplot(counts %>% filter(str_detect(Metric, "^Cigar")) %>% filter(!str_detect(Metric, "bp")), aes(x=ID, y=Value, fill=Squashed)) +
    geom_bar(stat="identity", position="dodge") +
    facet_grid(Metric~., scales="free") +
    theme_bw() +
    theme(axis.text.x = element_text(angle=45, hjust=1, size=7)) +
    ggtitle("Counts of variants detected in cigar string")
ggsave("plots/cigar_counts.png")

ggplot(counts %>% filter(str_detect(Metric, "^Split")), aes(x=ID, y=Value, fill=Squashed)) +
    geom_bar(stat="identity", position="dodge") +
    facet_grid(Metric~., scales="free") +
    theme_bw() +
    theme(axis.text.x = element_text(angle=45, hjust=1, size=7)) +
    ggtitle("Counts of variants detected from split alignments")
ggsave("plots/split_counts.png")

ggplot(counts %>% filter(str_detect(Metric, "^total")) %>% filter(!str_detect(Metric, "our")), aes(x=ID, y=Value, fill=Squashed)) +
    geom_bar(stat="identity", position="dodge") +
    facet_wrap("Metric") +
    theme_bw() +
    theme(axis.text.x = element_text(angle=45, hjust=1, size=7)) +
    ggtitle("Total 50bp+ SV calls combined cigar and split alignments")
ggsave("plots/total_sv_counts.png")

small_variants <- data.frame()
for (in_file in small_variants_files$file) {
    assembly <- gsub(".happy.extended.csv", "", basename(in_file))
    df <- read.table(file=in_file, header=TRUE, sep=",")
    df$Assembly <- assembly
    small_variants <- rbind(small_variants, df)
}
small_variants <- small_variants %>% filter(Subset == "*" & Subtype=="*") %>% filter(Filter == "PASS" | Filter=="*")
small_variants$Subset <- recode(small_variants$Subset, "*"="ALL")
small_variants$Type <- factor(small_variants$Type, levels=c("SNP", "INDEL"))
small_variants <- merge(small_variants, assemblies, by.x="Assembly", by.y="Name")
ggplot(small_variants, aes(x=ID, y=METRIC.Recall, fill=Squashed)) +
    geom_bar(stat="identity", position="dodge") +
    facet_grid(Type~.) +
    theme_bw() +
    theme(axis.text.x = element_text(angle=45, hjust=1, size=7)) +
    ggtitle("Hap.py sensitivity to truthset")
ggsave("plots/happy.sensitivity.png")

indels <- data.frame()
for (in_file in indels_files$file) {
    assembly <- gsub(".indels.indel.counts.horizontal.txt", "", basename(in_file))
    df <- read.table(file=in_file, header=FALSE, sep="\t", col.names=c("Metric", "Value"))
    df$Assembly <- assembly
    indels <- rbind(indels, df)
}
indels$slop <- 50
indels$slop[grepl("10 bp", indels$Metric)] <- 10
indels$slop[grepl("1 bp", indels$Metric)] <- 1
truth_name <- "GiaB_HG002"
indels$Metric <- recode(indels$Metric,
        "total GiaB_HG002 het indels"="totalTruthHet",
        "pairtopair for GiaB_HG002 indels vs our calls, 10 bp slop"="truthAllMatched",
        "pairtopair for GiaB_HG002 indels vs our calls, 50 bp slop"="truthAllMatched",
        "pairtopair for GiaB_HG002 indels vs our calls, 1 bp slop"="truthAllMatched",
        "pairtopair for GiaB_HG002 het to homalt, 50 bp slop"="truthHtToHomalt",
        "pairtopair for GiaB_HG002 het to homalt, 10 bp slop"="truthHetToHomalt",
        "pairtopair for GiaB_HG002 het to homalt, 1 bp slop"="truthHetToHomalt",
        "pairtopair for GiaB_HG002 het indels vs our calls, 50 bp slop"="truthHetMatched",
        "pairtopair for GiaB_HG002 het indels vs our calls, 10 bp slop"="truthHetMatched",
        "pairtopair for GiaB_HG002 het indels vs our calls, 1 bp slop"="truthHetMatched",
        "pairtopair for GiaB_HG002 het to het, 50 bp slop"="truthHetToHet",
        "pairtopair for GiaB_HG002 het to het, 10 bp slop"="truthHetToHet",
        "pairtopair for GiaB_HG002 het to het, 1 bp slop"="truthHetToHet")

for (assembly in unique(indels$Assembly)) {
    row <- indels[indels$Assembly==assembly & indels$Metric=="totalTruthHet",]
    row$slop <- 10
    indels <- rbind(indels, row)
    row$slop <- 1
    indels <- rbind(indels, row)
}
indels2 <- dcast(indels, Assembly+slop~Metric, value.var="Value")
indels2$truthHetToHomalt <- indels2$truthHetMatched - indels2$truthHetToHet
indels2$truthHetToHomref <- indels2$totalTruthHet - indels2$truthHetMatched
indels <- melt(indels2, id.vars=c("Assembly", "slop"), variable.name="Metric", value.name="Value")

indels <- merge(indels, assemblies, by.x="Assembly", by.y="Name")
ggplot(indels %>% filter(Metric == "truthHetToHomalt" | Metric == "truthHetToHet" | Metric == "truthHetToHomref"), aes(x=ID, y=Value, fill=Metric)) +
    geom_bar(stat="identity", position="stack") +
    facet_grid(slop~Squashed, scales="free") +
    theme_bw() +
    theme(axis.text.x = element_text(angle=45, hjust=1, size=5)) +
    ggtitle("Fates of heterozygous indels from truthset")
ggsave("plots/truthHetIndelFates.png", width=8, height=6, units="in")

ggplot(indels %>% filter(Metric == "truthHetToHomalt" | Metric == "truthHetToHomref"), aes(x=ID, y=Value, fill=Metric)) +
    geom_bar(stat="identity", position="stack") +
    facet_grid(slop~Squashed, scales="free") +
    theme_bw() +
    theme(axis.text.x = element_text(angle=45, hjust=1, size=6)) +
    ggtitle("Fates of heterozygous indels from truthset")
ggsave("plots/truthHetIndelFatesErrorsOnly.png")

sv <- data.frame()
for (in_file in sv_files$file) {
    assembly <- gsub(".sv.counts.horizontal.txt", "", basename(in_file))
    df <- read.table(file=in_file, header=FALSE, sep="\t", col.names=c("Metric", "Value"))
    df$Assembly <- assembly
    sv <- rbind(sv, df)
}
sv$Region <- "ALL"
sv$Region[grepl("str", sv$Metric)] <- "STR"
sv$Region[grepl("segDup", sv$Metric)] <- "SEGDUP"
sv$Region[grepl("nonRep", sv$Metric)] <- "NONREP"
sv$Metric <- gsub("str ", "", sv$Metric)
sv$Metric <- gsub("nonRep ", "", sv$Metric)
sv$Metric <- gsub("segDup ", "", sv$Metric)
sv$Metric <- gsub("str", "", sv$Metric)
sv$Metric <- gsub("nonRep", "", sv$Metric)
sv$Metric <- gsub("segDup", "", sv$Metric)

sv$Filter <- "ALL"
sv$Filter[grepl("PASS", sv$Metric)] <- "PASS"
sv$Metric <- gsub("PASS ", "", sv$Metric)
sv$Metric <- gsub("PASS", "", sv$Metric)

sv$Metric <- gsub("for ", "", sv$Metric)
sv$Metric <- gsub(", 50 bp slop", "", sv$Metric)
sv$Metric <- gsub(" $", "", sv$Metric)

sv$Overlap <- "pairtopair 50 bp slop"
sv$Overlap[grepl("50", sv$Metric)] <- "reciprocal 50% overlap"
sv$Overlap[grepl("10", sv$Metric)] <- "reciprocal 10% overlap"
sv$Metric <- gsub("%", "", sv$Metric)
sv$Metric <- gsub("50 ", "", sv$Metric)
sv$Metric <- gsub("10 ", "", sv$Metric)

sv$Overlap[sv$Metric=="total GiaB_HG002 bed"] <- "reciprocal 50% overlap"
rows <- sv[sv$Metric=="total GiaB_HG002 bed",]
rows$Overlap <- "reciprocal 10% overlap"
sv <- rbind(sv, rows)

sv2 <- dcast(sv, Assembly+Region+Filter+Overlap~Metric, value.var="Value")
sv2$sensitivity <- sv2$"pairtopair GiaB_HG002 vs combined calls"/sv2$"total GiaB_HG002 calls"
sv2$sensitivity[is.na(sv2$sensitivity)] <- sv2$"Reciprocal overlap GiaB_HG002"[is.na(sv2$sensitivity)]/sv2$"total GiaB_HG002 bed"[is.na(sv2$sensitivity)]
sv2$count <- sv2$"pairtopair GiaB_HG002 vs combined calls"
sv2$count[is.na(sv2$count)] <- sv2$"Reciprocal overlap GiaB_HG002"[is.na(sv2$count)]

sv2 <- merge(sv2, assemblies, by.x="Assembly", by.y="Name")
ggplot(sv2 %>% filter(Region != "ALL") %>% filter(Filter!= "ALL"), aes(x=ID, y=count, fill=Squashed)) +
    geom_bar(stat="identity", position="dodge") +
    facet_grid(Region~Overlap, scales="free") +
    theme_bw() +
    theme(axis.text.x = element_text(angle=45, hjust=1, size=6)) +
    ggtitle("SV sensitivity counts")
ggsave("plots/sv.sensitivityCounts.png", width=10, height=10, units="in")

ggplot(sv2 %>% filter(Region != "ALL") %>% filter(Filter != "ALL"), aes(x=ID, y=sensitivity, fill=Squashed)) +
    geom_bar(stat="identity", position="dodge") +
    facet_grid(Region~Overlap, scales="free") +
    theme_bw() +
    theme(axis.text.x = element_text(angle=45, hjust=1, size=6)) +
    ggtitle("SV sensitivity")
ggsave("plots/sv.sensitivity.png", width=10, height=10, units="in")

genome_cov <- data.frame()
cov_summary <- data.frame()
for (in_file in genome_cov_files$file) {
    assembly <- gsub(".genomecov.noGaps.hist", "", basename(in_file))
    cov <- read.table(file=in_file, header=FALSE, sep="\t", col.names=c("chrom", "cov", "bp"))
    cov <- cov %>% filter(chrom %in% c("chr1", "chr2", "chr3", "chr4", "chr5", "chr6", "chr7", "chr8", "chr9", "chr10", "chr11", "chr12", "chr13", "chr14", "chr15", "chr16", "chr17", "chr18", "chr19", "chr20", "chr21", "chr22", "chrX", "chrY"))
    gt9 <- cov %>% filter(cov > 9)
    lt10 <- cov %>% filter(cov < 10)
    gt9_zoomed <- gt9 %>% group_by(chrom) %>% summarize(bp=sum(bp))
    gt9_zoomed$cov <- 10
    cov_zoomed <- rbind(lt10, gt9_zoomed)
    cov_zoomed$Assembly <- assembly
    genome_cov <- rbind(genome_cov, cov_zoomed)
    ggplot(cov_zoomed, aes(x=cov, y=bp)) +
        geom_bar(stat="identity", position="dodge") +
        facet_wrap("chrom", scales="free") +
        scale_x_discrete(name="Coverage") +
        theme_bw() +
        ggtitle(assembly)
    ggsave(paste("plots", paste(assembly, "no_gaps_coverage.png", sep="."), sep="/"))
    covDip <- cov %>% filter(!grepl("X", chrom)) %>% filter(!grepl("Y", chrom))
    gt9 <- covDip %>% filter(cov > 9)
    lt10 <- covDip %>% filter(cov < 10)
    gt9_zoomed <- gt9 %>% group_by(chrom) %>% summarize(bp=sum(bp))
    gt9_zoomed$cov <- 10
    covDip_zoomed <- rbind(lt10, gt9_zoomed)
    summary_bp <- covDip_zoomed %>% group_by(cov) %>% summarize(sum=as.numeric(sum(as.numeric(bp))))
    total_bp <- sum(summary_bp$sum)
    summary_bp$percent <- summary_bp$sum/total_bp
    summary_bp$Assembly <- assembly
    cov_summary <- rbind(cov_summary, summary_bp)

    ggplot(summary_bp, aes(x=cov, y=percent)) +
        geom_bar(stat="identity", position="dodge") +
        scale_x_discrete(name="Coverage") +
        theme_bw() +
        ggtitle(assembly)
    ggsave(paste("plots", paste(assembly, "no_gaps_coverage_summary.png", sep="."), sep="/"))
}
genome_cov <- merge(genome_cov, assemblies, by.x="Assembly", by.y="Name")
ggplot(genome_cov, aes(x=cov, y=bp)) +
    geom_boxplot() +
    facet_grid(Squashed~chrom, scales="free") +
    scale_x_discrete(name="Coverage") +
    theme_bw() +
    ggtitle("Coverage by chromosome, excludes ref gaps")
ggsave(paste("plots", "no_gaps_coverage.png", sep="/"))

cov_summary <- merge(cov_summary, assemblies, by.x="Assembly", by.y="Name")
cov_summary$cov <- as.character(cov_summary$cov)
cov_summary$cov[cov_summary$cov=="10"] <- "10+"
cov_summary$cov <- factor(cov_summary$cov, levels=c("0", "1", "2", "3", "4", "5", "6", "7", "8", "9", "10+"))
ggplot(cov_summary, aes(x=cov, y=percent)) +
    geom_boxplot() +
    facet_grid(.~Squashed, scales="free") +
    scale_x_discrete(name="Coverage") +
    theme_bw() +
    ggtitle("Coverage for chr1-22, excludes ref gaps")
ggsave(paste("plots", "no_gaps_coverage_summary.png", sep="/"))

for (in_file in nonRep_lengths_files$file) {
    assembly <- gsub(".ours.nonRep.lengths.txt", "", basename(in_file))
    sv_nonrep <- read.table(in_file, header=FALSE, col.names=c("type", "length"))
    sv_nonrep$region <- "nonRep"
    sv_str <- read.table(gsub("nonRep", "str", in_file), header=FALSE, col.names=c("type", "length"))
    sv_str$region <- "str"
    sv_segdup <- read.table(gsub("nonRep", "segDup", in_file), header=FALSE, col.names=c("type", "length"))
    sv_segdup$region <- "segDup"
    sv <- rbind(sv_nonrep, sv_str, sv_segdup) %>% filter(type!="GAP")

    sv$length_zoom <- sv$length
    sv$length_zoom[sv$length_zoom < -10000000] <- -1050.5
    sv$length_zoom[sv$length_zoom < -1000000] <- -1040.5
    sv$length_zoom[sv$length_zoom < -100000] <- -1030.5
    sv$length_zoom[sv$length_zoom < -10000] <- -1020.5
    sv$length_zoom[sv$length_zoom < -1000 & sv$length_zoom != -1050.5 & sv$length_zoom != -1040.5 & sv$length_zoom != -1030.5 & sv$length_zoom != -1020.5] <- -1010.5
    sv$length_zoom[sv$length_zoom > 1000 & sv$length_zoom < 10000] <- 1010.5
    sv$length_zoom[sv$length_zoom > 10000 & sv$length_zoom < 100000] <- 1020.5
    sv$length_zoom[sv$length_zoom > 100000 & sv$length_zoom < 1000000] <- 1030.5
    sv$length_zoom[sv$length_zoom > 1000000 & sv$length_zoom < 10000000] <- 1040.5
    sv$length_zoom[sv$length_zoom > 10000000] <- 1050.5
    sv$length_zoom[sv$length_zoom==-1050.5] <- -1050
    sv$length_zoom[sv$length_zoom==-1040.5] <- -1040
    sv$length_zoom[sv$length_zoom==-1030.5] <- -1030
    sv$length_zoom[sv$length_zoom==-1020.5] <- -1020
    sv$length_zoom[sv$length_zoom==-1010.5] <- -1010
    sv$length_zoom[sv$length_zoom==1010.5] <- 1010
    sv$length_zoom[sv$length_zoom==1020.5] <- 1020
    sv$length_zoom[sv$length_zoom==1030.5] <- 1030
    sv$length_zoom[sv$length_zoom==1040.5] <- 1040
    sv$length_zoom[sv$length_zoom==1050.5] <- 1050
    ggplot(sv, aes(x=length_zoom)) +
        geom_histogram(binwidth=10) +
        facet_grid(region~type, scales="free") +
        theme_bw() +
        theme(axis.text.x = element_text(angle = 45, hjust = 1, size=7)) +
        ggtitle(assembly)
    ggsave(file=paste("plots", paste(assembly, "sv_lengths.png", sep="."), sep="/"))
}

svCounts <- data.frame()
for (in_file in nonRep_types_files$file) {
    assembly <- gsub(".ours.nonRep.typeCount.txt", "", basename(in_file))
    for (type in c("nonRep", "str", "segDup")) {
        df <- read.table(gsub("nonRep", type, in_file), header=FALSE, sep=" ", col.names=c("Value", "Metric"))
        df$Assembly <- assembly
        df$Type <- type
        svCounts <- rbind(svCounts, df)
    }
}
svCounts$Sample <- gsub("\\..*", "", svCounts$Assembly)
svCounts <- merge(svCounts, populations, by="Sample")
svCounts <- merge(svCounts, assemblies, by.x="Assembly", by.y="Name")

ggplot(svCounts %>% filter(Type=="nonRep") %>% filter(!str_detect(Metric, "^G")), aes(x=ID, y=Value, fill=Superpopulation)) +
    geom_bar(stat="identity", position="dodge") +
    facet_wrap("Metric", scales="free") +
    theme_bw() +
    theme(axis.text.x = element_text(angle=45, hjust=1, size=6)) +
    ggtitle("Non-repetitive SV counts by type")
    ggsave("plots/nonrepSvCountsPopulation.png", width=12, height=8, units="in")

ggplot(svCounts %>% filter(Type=="nonRep") %>% filter(!str_detect(Metric, "^G")), aes(x=ID, y=Value, fill=Squashed)) +
    geom_bar(stat="identity", position="dodge") +
    facet_wrap("Metric", scales="free") +
    theme_bw() +
    theme(axis.text.x = element_text(angle=45, hjust=1, size=6)) +
    ggtitle("Non-repetitive SV counts by type")
    ggsave("plots/nonrepSvCountsSquashed.png", width=12, height=8, units="in")

ggplot(svCounts %>% filter(Type=="str") %>% filter(!str_detect(Metric, "^G")), aes(x=ID, y=Value, fill=Superpopulation)) +
    geom_bar(stat="identity", position="dodge") +
    facet_wrap("Metric", scales="free") +
    theme_bw() +
    theme(axis.text.x = element_text(angle=45, hjust=1, size=6)) +
    ggtitle("STR SV counts by type")
ggsave("plots/strSvCountsPopulation.png", width=12, height=8, units="in")

ggplot(svCounts %>% filter(Type=="str") %>% filter(!str_detect(Metric, "^G")), aes(x=ID, y=Value, fill=Squashed)) +
    geom_bar(stat="identity", position="dodge") +
    facet_wrap("Metric", scales="free") +
    theme_bw() +
    theme(axis.text.x = element_text(angle=45, hjust=1, size=6)) +
    ggtitle("STR SV counts by type")
ggsave("plots/strSvCountsSquashed.png", width=12, height=8, units="in")

ggplot(svCounts %>% filter(Type=="segDup") %>% filter(!str_detect(Metric, "^G")), aes(x=ID, y=Value, fill=Superpopulation)) +
    geom_bar(stat="identity", position="dodge") +
    facet_wrap("Metric", scales="free") +
    theme_bw() +
    theme(axis.text.x = element_text(angle=45, hjust=1, size=6)) +
    ggtitle("SegDup SV counts by type")
ggsave("plots/segDupSvCountsPopulation.png", width=12, height=8, units="in")

ggplot(svCounts %>% filter(Type=="segDup") %>% filter(!str_detect(Metric, "^G")), aes(x=ID, y=Value, fill=Squashed)) +
    geom_bar(stat="identity", position="dodge") +
    facet_wrap("Metric", scales="free") +
    theme_bw() +
    theme(axis.text.x = element_text(angle=45, hjust=1, size=6)) +
    ggtitle("SegDup SV counts by type")
ggsave("plots/segDupSvCountsSquashed.png", width=12, height=8, units="in")

het <- data.frame()
for (in_file in het_fates_files$file) {
    assembly <- gsub(".happy.het.counts.horizontal.txt", "", basename(in_file))
    df <- read.table(in_file, header=FALSE, col.names=c("Type", "Origin", "Fate", "Value"))
    df$Assembly <- assembly
    het <- rbind(het, df)
}
het$Type <- factor(het$Type, levels=c("SNP", "INDEL"))
het <- merge(het, assemblies, by.x="Assembly", by.y="Name")
het$ID[het$Squashed] <- paste(het$ID[het$Squashed], "*")
ggplot(het, aes(x=ID, y=Value, fill=Fate)) +
    geom_bar(stat="identity", position="stack") +
    facet_grid(Type~Origin, scales="free") +
    theme_bw() +
    theme(axis.text.x = element_text(angle=45, hjust=1, size=7)) +
    ggtitle("Fate of variants from truthset")
ggsave("plots/het_fates.png")

rep <- data.frame()
for (in_file in str_venn_files$file) {
    assembly <- gsub(".sv.str.VennInput.txt", "", basename(in_file))
    df <- read.table(in_file, header=FALSE, row.names=c("A", "B", "C", "AandB", "BandC", "AandC", "AandBandC"))
    df2 <- transpose(df)
    colnames(df2) <- rownames(df)
    df2$Assembly <- assembly
    df2$Region <- "str"
    rep <- rbind(rep, df2)
    df <- read.table(gsub("str", "segDup", in_file), header=FALSE, row.names=c("A", "B", "C", "AandB", "BandC", "AandC", "AandBandC"))
    df2 <- transpose(df)
    colnames(df2) <- rownames(df)
    df2$Assembly <- assembly
    df2$Region <- "segDup"
    rep <- rbind(rep, df2)
}
rep$totalGiaBPass <- rep$A
rep$totalGiaB <- rep$B
rep$totalOurs <- rep$C
rep$sensitivity <- rep$BandC/rep$totalGiaB
rep$sensitivityPass <- rep$AandC/rep$totalGiaBPass
rep$count <- rep$BandC
rep$countPass <- rep$AandC
rep2 <- melt(rep, id.vars=c("Assembly", "Region"))
rep2 <- merge(rep2, assemblies, by.x="Assembly", by.y="Name")
ggplot(rep2 %>% filter(variable=="sensitivityPass"), aes(x=ID, y=value, fill=Squashed)) +
    geom_bar(stat="identity", position="dodge") +
    facet_grid(Region~., scales="free") +
    theme_bw() +
    theme(axis.text.x = element_text(angle=45, hjust=1, size=7)) +
    ggtitle("Sensitivity to repetitive element detection")
ggsave("plots/truthset.repetitive.sensitivity.png")

rep$unfound <- rep$B - rep$BandC
rep$unfoundPass <- rep$A - rep$AandC
rep2 <- melt(rep, id.vars=c("Assembly", "Region"))
rep2$found <- FALSE
rep2$found[rep2$variable=="count" | rep2$variable=="countPass"] <- TRUE
rep2$plot <- FALSE
rep2$plot[rep2$variable=="countPass" | rep2$variable=="unfoundPass"] <- TRUE
rep2$Filter <- "PASS"
ggplot(rep2 %>% filter(rep2$plot==TRUE), aes(x=Assembly, y=value, fill=found)) +
    geom_bar(stat="identity", position="stack") +
    facet_grid(Region~., scales="free") +
    theme_bw() +
    theme(axis.text.x = element_text(angle=45, hjust=1, size=7)) +
    ggtitle("Repetitive element detection counts")
ggsave("plots/truthset.repetitive.count.png")
