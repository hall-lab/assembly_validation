version 1.0
import "analyze_assembly.wdl" as analyze

workflow PlotAssemblies {
    input {
        File populations
        File assembly_list
        File ref
        File ref_index
	File gap_file
	File str_track
	File seg_dup_track
        File truth_list_snp_indel
        File truth_list_confident_region
        File truth_list_sv_vcf
    }
    Array[Array[File]] assemblies = read_tsv(assembly_list)
    Map[String, File] truthsets_snp_indel = read_map(truth_list_snp_indel)
    Map[String, File] truthsets_confident_region = read_map(truth_list_confident_region)
    Map[String, File] truthsets_sv_vcf = read_map(truth_list_sv_vcf)

    scatter (assembly in assemblies) {
	String truth_name = assembly[3]
	call analyze.AnalyzeAssembly as analyze_assembly {
	    input:
		contigs1=assembly[1],
		contigs2=assembly[2],
		ref=ref,
		ref_index=ref_index,
		truth_vcf_snp_indel=if truth_name != "NONE" then truthsets_snp_indel[truth_name] else "None",
		truth_confident_region=if truth_name != "NONE" then truthsets_confident_region[truth_name] else "None",
		truth_vcf_sv=if truth_name != "NONE" then truthsets_sv_vcf[truth_name] else "None",
		gap_file=gap_file,
		str_track=str_track,
		seg_dup_track=seg_dup_track,
		output_prefix=assembly[0],
		sample_name=assembly[0],
		truth_name=truth_name
	}
	
    }

    call make_plots {
        input:
            populations=populations,
            small_variants=select_all(analyze_assembly.small_variants),
            indels=select_all(analyze_assembly.indels),
            sv=select_all(analyze_assembly.sv),
            genome_cov=analyze_assembly.genome_cov,
            nonRep_lengths=analyze_assembly.nonRep_lengths,
            str_lengths=analyze_assembly.str_lengths,
            segDup_lengths=analyze_assembly.segDup_lengths,
            counts=analyze_assembly.counts,
            nonRep_types=analyze_assembly.nonRep_types,
            str_types=analyze_assembly.str_types,
            segDup_types=analyze_assembly.segDup_types,
            het_fates=select_all(analyze_assembly.het_fates),
            str_venn=select_all(analyze_assembly.str_venn),
            segDup_venn=select_all(analyze_assembly.segDup_venn)
    }

    output {
    	File plot_tarball = make_plots.plots_tar
    }
}

task make_plots {
    input {
        Array[File] small_variants
        Array[File] indels
        Array[File] sv
        Array[File] genome_cov
        Array[File] nonRep_lengths
        Array[File] str_lengths
        Array[File] segDup_lengths
        Array[File] counts
        Array[File] nonRep_types
        Array[File] str_types
        Array[File] segDup_types
        Array[File] het_fates
        Array[File] str_venn
        Array[File] segDup_venn
        File populations
        File small_variants_fof = write_lines(small_variants)
        File indels_fof = write_lines(indels)
        File sv_fof = write_lines(sv)
        File genome_cov_fof = write_lines(genome_cov)
        File nonRep_lengths_fof = write_lines(nonRep_lengths)
        File str_lengths_fof = write_lines(str_lengths)
        File segDup_lengths_fof = write_lines(segDup_lengths)
        File counts_fof = write_lines(counts)
        File nonRep_types_fof = write_lines(nonRep_types)
        File str_types_fof = write_lines(str_types)
        File segDup_types_fof = write_lines(segDup_types)
        File het_fates_fof = write_lines(het_fates)
        File str_venn_fof = write_lines(str_venn)
        File segDup_venn_fof = write_lines(segDup_venn)
    }
    command <<<
        set -exo pipefail
        mkdir -p plots
        RSCRIPT=/usr/local/bin/Rscript
        ${RSCRIPT} /opt/hall-lab/make_plots.R ~{populations} ~{small_variants_fof} ~{indels_fof} ~{sv_fof} ~{genome_cov_fof} ~{nonRep_lengths_fof} ~{str_lengths_fof} ~{segDup_lengths_fof} ~{counts_fof} ~{nonRep_types_fof} ~{str_types_fof} ~{segDup_types_fof} ~{het_fates_fof} ~{str_venn_fof} ~{segDup_venn_fof}
        tar -czf plots.tar.gz plots
    >>>
    runtime {
        docker: "apregier/plot_assemblies@sha256:32e0f400595dc98c05a2fa31023e9558bd51c71e60e2069b48222b561d459668"
        memory: "4 GB"
    }
    output {
        File plots_tar="plots.tar.gz"
    }
}
