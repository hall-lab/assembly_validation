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
        File truth_list_sv_confident_region
        File jupyter_notebook_template
    }
    Array[Array[File]] assemblies = read_tsv(assembly_list)
    Map[String, File] truthsets_snp_indel = read_map(truth_list_snp_indel)
    Map[String, File] truthsets_confident_region = read_map(truth_list_confident_region)
    Map[String, File] truthsets_sv_vcf = read_map(truth_list_sv_vcf)
    Map[String, File] truthsets_sv_confident_region = read_map(truth_list_sv_confident_region)

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
        truth_confident_region_sv=if truth_name != "NONE" then truthsets_sv_confident_region[truth_name] else "None",
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
            segDup_venn=select_all(analyze_assembly.segDup_venn),
            cigar_indel_lengths=analyze_assembly.cigar_indel_lengths,
            assembly_file=assembly_list, 
            template=jupyter_notebook_template
    }

    output {
        File jupyter_notebook = make_plots.jupyter_notebook
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
        Array[File] cigar_indel_lengths
        File populations
        File assembly_file
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
        File cigar_indel_lengths_fof = write_lines(cigar_indel_lengths)
    }
    command <<<
        set -exo pipefail
        mkdir -p plots
        SCRIPT=/opt/hall-lab/customize_notebook.py
        echo "counts_horizontal	~{counts_fof}" > notebook_inputs.tsv
        echo "str_typeCount	~{str_types_fof}" >> notebook_inputs.tsv
        echo "nonRep_typeCount	~{nonRep_types_fof}" >> notebook_inputs.tsv
        echo "nonRep_lengths	~{nonRep_lengths_fof}" >> notebook_inputs.tsv
        echo "sv_str_VennInput	~{str_venn_fof}" >> notebook_inputs.tsv
        echo "str_lengths	~{str_lengths_fof}" >> notebook_inputs.tsv
        echo "happy_extended	~{small_variants_fof}" >> notebook_inputs.tsv
        echo "genomecov_noGaps	~{genome_cov_fof}" >> notebook_inputs.tsv
        echo "segDup_typeCount	~{segDup_types_fof}" >> notebook_inputs.tsv
        echo "indel_lengths_horizontal	~{cigar_indel_lengths_fof}" >> notebook_inputs.tsv
        echo "happy_het_counts_horizontal	~{het_fates_fof}" >> notebook_inputs.tsv
        echo "segDup_lengths	~{segDup_lengths_fof}" >> notebook_inputs.tsv
        echo "sv_segDup_VennInput	~{segDup_venn_fof}" >> notebook_inputs.tsv
        echo "indel_counts	~{indels_fof}" >> notebook_inputs.tsv
        echo "sv_counts_horizontal	~{sv_fof}" >> notebook_inputs.tsv
        $SCRIPT --tsv notebook_inputs.tsv --template ~{template} --output new.ipynb
    >>>
    runtime {
        docker: "apregier/plot_assemblies@sha256:436cb951a37b4fb953d98d9fdf1aad556db538f7df4b17a04ea2354d55ea94e1"
        memory: "4 GB"
    }
    output {
        File jupyter_notebook="new.ipynb"
    }
}
