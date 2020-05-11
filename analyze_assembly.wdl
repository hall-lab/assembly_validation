version 1.0

workflow AnalyzeAssembly {
    input {
        File contigs1
        File contigs2
        File ref
        File ref_index
        File? truth_vcf_snp_indel
        File? truth_confident_region
        File? truth_vcf_sv
        File? truth_confident_region_sv
        File gap_file
        File str_track
        File seg_dup_track
        String output_prefix
        String sample_name
        String truth_name
    }

    call align_contigs {
        input:
            contigs1=contigs1,
            contigs2=contigs2,
            ref=ref,
            output_prefix=output_prefix,
            sample_name=sample_name,
            gap_file=gap_file,
            str_track=str_track,
            seg_dup_track=seg_dup_track
    }
    call get_indel_lengths {
        input:
            vcf=align_contigs.genotyped_vcf,
            output_prefix=output_prefix
    }
    if (truth_name != "NONE") {
        call compare_indels {
            input:
                vcf=align_contigs.genotyped_vcf,
                truth=truth_vcf_snp_indel,
                output_prefix=output_prefix,
                sample_name=sample_name,
                truth_name=truth_name,
                str_bed=align_contigs.str_bed,
                seg_dup_bed=align_contigs.seg_dup_bed
        }
        call compare_svs {
            input:
                bedpe=align_contigs.bedpe,
                truth=truth_vcf_sv,
                truth_confident_region=truth_confident_region_sv,
                output_prefix=output_prefix,
                sample_name=sample_name,
                truth_name=truth_name,
                str_bed=align_contigs.str_bed,
                seg_dup_bed=align_contigs.seg_dup_bed
        }
        call compare_small {
            input:
                ref=ref,
                ref_index=ref_index,
                truth=truth_vcf_snp_indel,
                vcf=align_contigs.genotyped_vcf,
		vcf_index=align_contigs.genotyped_vcf_index,
                confident_region=truth_confident_region,
                output_prefix=output_prefix
        }
        call calculate_het_fates {
            input:
                vcf=compare_small.happy_vcf,
                output_prefix=output_prefix
        }
    }

    output {
            File? small_variants=compare_small.counts
            File? indels=compare_indels.counts
            File? sv=compare_svs.counts
            File genome_cov=align_contigs.genome_cov
            File nonRep_lengths=align_contigs.nonRep_lengths
            File str_lengths=align_contigs.str_lengths
            File segDup_lengths=align_contigs.segDup_lengths
            File counts=align_contigs.counts
            File nonRep_types=align_contigs.nonRep_types
            File str_types=align_contigs.str_types
            File segDup_types=align_contigs.segDup_types
            File cigar_indel_lengths=get_indel_lengths.cigar_indel_lengths
            File? het_fates=calculate_het_fates.het_counts
            File? str_venn=compare_svs.str_venn
            File? segDup_venn=compare_svs.segDup_venn
    }
}

task compare_indels {
    input {
        File vcf
        File? truth
        File str_bed
        File seg_dup_bed
        String output_prefix
        String sample_name
        String truth_name
    }
    command <<<
        set -exo pipefail
        /bin/bash /opt/hall-lab/scripts/compareIndels.sh ~{truth} ~{vcf} ~{output_prefix}.indels ~{sample_name} ~{truth_name} ~{seg_dup_bed} ~{str_bed}
    >>>
    runtime {
        memory: "16G"
        docker: "apregier/analyze_assemblies@sha256:edf94bd952180acb26423e9b0e583a8b00d658ac533634d59b32523cbd2a602a"
    }
    output {
        File counts="${output_prefix}.indels.indel.counts.horizontal.txt"
    }
}

task compare_svs {
    input {
        File bedpe
        File? truth
        File? truth_confident_region
        File str_bed
        File seg_dup_bed
        String output_prefix
        String sample_name
        String truth_name
    }
    command <<<
        set -exo pipefail
        /bin/bash /opt/hall-lab/scripts/compareSVs.sh ~{truth} ~{truth_confident_region} ~{bedpe} ~{output_prefix}.sv ~{sample_name} ~{truth_name} ~{seg_dup_bed} ~{str_bed}
    >>>
    runtime {
        memory: "8G"
        docker: "apregier/analyze_assemblies@sha256:edf94bd952180acb26423e9b0e583a8b00d658ac533634d59b32523cbd2a602a"
    }
    output {
        File counts="${output_prefix}.sv.counts.horizontal.txt"
        File str_venn="${output_prefix}.sv.str.VennInput.txt"
        File segDup_venn="${output_prefix}.sv.segDup.VennInput.txt"
    }
}

task compare_small {
    input {
        File ref
        File ref_index
        File? truth
        File vcf
	File vcf_index
        File? confident_region
        String output_prefix
    }
    command <<<
        HGREF=~{ref} /opt/hap.py/bin/hap.py ~{truth} ~{vcf} -o ~{output_prefix}.happy \
        -f ~{confident_region} \
        --location chr1,chr2,chr3,chr4,chr5,chr6,chr7,chr8,chr9,chr10,chr11,chr12,chr13,chr14,chr15,chr16,chr17,chr18,chr19,chr20,chr21,chr22 \
        --threads 4 \
        --preprocess-truth
    >>>
    runtime {
        memory: "32G"
        docker: "pkrusche/hap.py"
    }
    output {
        File counts="${output_prefix}.happy.extended.csv"
        File happy_vcf="${output_prefix}.happy.vcf.gz"
    }
}

task calculate_het_fates {
    input {
        File vcf
        String output_prefix
    }
    command <<<
        /bin/bash /opt/hall-lab/scripts/calculate_het_fates.sh ~{vcf} ~{output_prefix}.happy
    >>>
    runtime {
        memory: "4G"
        docker: "apregier/analyze_assemblies@sha256:edf94bd952180acb26423e9b0e583a8b00d658ac533634d59b32523cbd2a602a"
    }
    output {
        File het_counts="${output_prefix}.happy.het.counts.horizontal.txt"
    }
}

task get_indel_lengths {
    input {
        File vcf
        String output_prefix
    }
    command <<<
        set -exo pipefail
        touch tmp.2
        /opt/hall-lab/python-2.7.15/bin/python /opt/hall-lab/scripts/sv_lengths.py -v ~{vcf} -o ~{output_prefix}.indel_lengths.horizontal.txt
    >>>
    runtime {
        memory: "4G"
        docker: "apregier/analyze_assemblies@sha256:edf94bd952180acb26423e9b0e583a8b00d658ac533634d59b32523cbd2a602a"
    }
    output {
        File cigar_indel_lengths="${output_prefix}.indel_lengths.horizontal.txt"
    }
}

task align_contigs {
    input {
            File contigs1
            File contigs2
            File ref
            File str_track
            File seg_dup_track
	    File gap_file
            String output_prefix
            String sample_name
    }
    command <<<
        set -exo pipefail
        /bin/bash /opt/hall-lab/scripts/align_contigs_to_ref.sh ~{contigs1} ~{contigs2} ~{ref} ~{output_prefix} ~{sample_name} ~{gap_file} ~{str_track} ~{seg_dup_track}
    >>>
    runtime {
        memory: "64G"
        docker: "apregier/analyze_assemblies@sha256:edf94bd952180acb26423e9b0e583a8b00d658ac533634d59b32523cbd2a602a"
    }
    output {
        File genotyped_vcf="${output_prefix}.loose.genotyped.vcf.gz"
        File genotyped_vcf_index="${output_prefix}.loose.genotyped.vcf.gz.tbi"
        File bedpe="${output_prefix}.ours.bedpe"
        File genome_cov="${output_prefix}.genomecov.noGaps.hist"
        File nonRep_lengths="${output_prefix}.ours.nonRep.lengths.txt"
        File str_lengths="${output_prefix}.ours.str.lengths.txt"
        File segDup_lengths="${output_prefix}.ours.segDup.lengths.txt"
        File counts="${output_prefix}.counts.horizontal.txt"
        File nonRep_types="${output_prefix}.ours.nonRep.typeCount.txt"
        File str_types="${output_prefix}.ours.str.typeCount.txt"
        File segDup_types="${output_prefix}.ours.segDup.typeCount.txt"
        File str_bed="${output_prefix}.str.merged.bed"
        File seg_dup_bed="${output_prefix}.segDup.bed"
    }
}
