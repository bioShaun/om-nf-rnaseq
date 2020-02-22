// TODO LIST
// - different libtype choice
// - optional analysis flag: rmats / snp / ...

def helpMessage() {
    log.info """

    A RNAseq analysis pipeline

    Usage:

    ==================================================================
    
    Mandatory arguments:
      --proj_name                   Project name for output name
      --reads                       Path to input data (must be surrounded with quotes)
      --sample_group                Sample vs Group file, tab seperated
      --outdir                      Path to analysis results 
      --db_name

    Optional:
      --rm_rrna                     Remove rRNA sequence in reads
      --qc                          Run qc analysis
      --pipeline                    Run whole rnaseq pipeline
      --kegg_abbr                   kegg abbr for kegg analysis

    """.stripIndent()
}

// workflow internal path&files
script_dir = file("${baseDir}/script/")
cfg_dir = file("${baseDir}/cfg/")

 // Show help emssage
params.help = false
if (params.help){
    helpMessage()
    exit 0
}


def check_ref_exist = {file_path, file_type ->
    if (file_path) {
        file_path = file(file_path)
        if( !file_path.exists() ) exit 1, "${file_type} file not found: ${file_path}"
        return file_path
    } else {
        exit 1, "${file_type} is required!"
    }
}

// annotations
params.genome = "${params.genomes_base}/${params.db_name}/${params.fasta_name}"
params.genome_fai = "${params.genomes_base}/${params.db_name}/${params.fasta_name}.fai"
params.gtf = "${params.genomes_base}/${params.db_name}/${params.gtf_name}"
params.gene_length = "${params.genomes_base}/${params.db_name}/${params.gene_len_name}"
params.split_gtf_dir = "${params.genomes_base}/${params.db_name}/${params.split_gtf_name}"
params.ref_flat = "${params.genomes_base}/${params.db_name}/${params.ref_flat_name}"
params.bed12 = "${params.genomes_base}/${params.db_name}/${params.bed12_name}"
params.chr_size = "${params.genomes_base}/${params.db_name}/${params.chr_size_name}"
params.chr_window = "${params.genomes_base}/${params.db_name}/${params.chr_window_name}.${params.window_size}.bed"  
params.star_index = "${params.genomes_base}/${params.db_name}/${params.star_index_name}"
params.kegg_anno = "${params.genomes_base}/${params.db_name}/${params.kegg_abbr}.${params.kegg_anno_name}"
params.go_anno = "${params.genomes_base}/${params.db_name}/${params.go_name}"
params.gene_ann = "${params.genomes_base}/${params.db_name}/${params.gene_des_name}"
params.venn = false


if (!params.kegg_backgroud) {
    kegg_backgroud = params.kegg_abbr
} else {
    kegg_backgroud = params.kegg_backgroud
}


if (params.pipeline) {
    sample_group = check_ref_exist(params.sample_group, 'Sample vs Group file')
} else {
    sample_group = ''
}

if (params.venn) {
    venn_list = file(params.venn)
} else {
    venn_list = file('')
}

// check ref
genome = check_ref_exist(params.genome, 'genome fasta')
genome_fai = check_ref_exist(params.genome_fai, 'genome fasta index')
gtf = check_ref_exist(params.gtf, 'gtf file')
bed12 = check_ref_exist(params.bed12, 'bed file')
ref_flat = check_ref_exist(params.ref_flat, 'ref flat file')
star_index = check_ref_exist(params.star_index, 'STAR index')
split_gtf_dir = check_ref_exist(params.split_gtf_dir, 'split gtf directory')
gene_ann = check_ref_exist(params.gene_ann, 'Gene annotation')
go_anno = check_ref_exist(params.go_anno, 'GO')
gene_length = check_ref_exist(params.gene_length, 'Gene Length')
kegg_anno = check_ref_exist(params.kegg_anno, 'KEGG')


// Prepare analysis fastq files
Channel
    .fromFilePairs("${params.reads}")
    .ifEmpty { exit 1, "Cannot find any reads matching: ${params.reads}\n!" }
    .set { raw_fq_files }


// contrast file

if (params.contrast) {
    contrast_file = file(params.contrast)
} else {
    process generate_contrast {
        publishDir "${params.outdir}/${params.proj_name}/configure" , mode: 'copy'

        input:
        file sample_group from sample_group
        
        output:
        file 'contrast.ini' into contrast_file
        
        when:
        params.pipeline

        script:
        """
        python ${script_dir}/utils/make_analysis_compare.py \\
            generate_contrast \\
            --sample-group ${sample_group} \\
            --contrast-file contrast.ini
        """
    }
}

//
if (params.rm_rrna) {
    rrna_idx = Channel
                    .fromPath("${params.rrna_index}.*.bt2")
                    .ifEmpty { exit 1, "rRNA index not found: ${params.rrna_index}" }
}
rrna_idx_base = file(params.rrna_index).getName()


// software version
process get_software_versions {

    output:
    file 'software_versions_mqc.yaml' into software_versions_yaml, software_versions_yaml_qc

    script:
    """
    echo $workflow.nextflow.version &> v_nextflow.txt
    fastp --version &> v_fastp.txt
    STAR --version &> v_star.txt
    samtools --version &> v_samtools.txt
    picard MarkDuplicates --version &> v_markduplicates.txt  || true
    multiqc --version &> v_multiqc.txt
    /public/software/stringtie/stringtie-2.0.6.Linux_x86_64/stringtie --version &> v_stringtie.txt
    gffcompare --version &> v_gffcompare.txt
    python ${script_dir}/utils/scrape_software_versions.py &> software_versions_mqc.yaml
    """
}

/*
 * Fastp
 */

process fastp {

    tag "${name}"

    publishDir "${params.outdir}/${params.proj_name}/result/data_info/detail/${name}", mode: 'copy',
        saveAs: {filename -> 
            if (filename.indexOf("fq.gz") > 0)  null
            else "${filename}" 
            }   

    cpus = 4

    input:
    set name, file(reads) from raw_fq_files

    output:
    file "*.fq.gz" into clean_reads, quant_reads
    file "${name}.json" into fastp_json_analysis, fastp_json_report
    file "${name}.html" into fastp_html

    script:
    """
    fastp \\
        --in1 ${reads[0]} \\
        --in2 ${reads[1]} \\
        --out1 ${name}.R1.clean.fq.gz \\
        --out2 ${name}.R2.clean.fq.gz \\
        --json ${name}.json \\
        --html ${name}.html    
    """
}  


/*
* fq quality control
*/
process reads_qc_summary {

    publishDir "${params.outdir}/${params.proj_name}/result/data_info/", mode: 'copy'

    input:
    file 'fastp_json/*' from fastp_json_analysis.collect()
    
    output:
    file "data_summary.csv" into data_summary, data_summary_qc
    file "reads_gc" into gc_plot, gc_plot_qc
    file "reads_quality" into rq_plot, rq_plot_qc
    file "reads_filter" into rf_plot, rf_plot_qc
    
    script:
    """
    python ${script_dir}/reads_qc/extract_fastp_info.py \\
        --fastp-dir fastp_json \\
        --outdir .

    /public/software/R/R-3.5.1/executable/bin/Rscript ${script_dir}/reads_qc/gc_plot.R \\
        --gc_dir reads_gc

    /public/software/R/R-3.5.1/executable/bin/Rscript ${script_dir}/reads_qc/reads_quality_plot.R \\
        --rq_dir reads_quality

    /public/software/R/R-3.5.1/executable/bin/Rscript ${script_dir}/reads_qc/filter_stats.R \\
        --filter_dir reads_filter   
    """
}

process star_mapping {
    tag "${sample_name}"

    publishDir "${params.outdir}/${params.proj_name}/result/alignment/${sample_name}", mode: 'copy',
        saveAs: {filename -> 
            if (filename.indexOf("Log.final.out") > 0)  "${sample_name}.star.log.txt"
            else if (filename.indexOf("mapping_rate") > 0) "${filename}"
            else null 
            }

    cpus = 40

    queue 'bm'

    input:
    file reads from clean_reads
    file index from star_index
    
    output:
    file "${sample_name}.bam" into unsort_bam
    file "${sample_name}.Log.final.out" into star_log, star_log_qc
    file "${sample_name}.mapping_rate*" into mapping_plt, mapping_plt_qc

    script:
    sample_name = reads[0].toString() - '.R1.clean.fq.gz'
    if (params.mrna_strand == 'unstranded') {
        strand_field = "--outSAMstrandField intronMotif"
    } else {
        strand_field = ""
    }
    """
    STAR    --genomeDir ${index} \\
            --readFilesIn ${reads[0]} ${reads[1]} \\
            --runThreadN ${task.cpus} \\
            --outFileNamePrefix ${sample_name}. \\
            --outSAMattrRGline  ID:${sample_name} CN:TCuni SM:${sample_name} LB:${sample_name} PI:350 PL:Illumina \\
            --readFilesCommand zcat \\
            --outSAMtype BAM Unsorted  \\
            --limitBAMsortRAM 100000000000 \\
            --outFilterType BySJout \\
            --outFilterMultimapNmax 20 \\
            --alignSJoverhangMin 8 \\
            --alignSJDBoverhangMin 1 \\
            --outFilterMismatchNmax 999 \\
            --alignIntronMin 20 \\
            --alignIntronMax 1000000 \\
            --alignMatesGapMax 1000000 \\
            --chimSegmentMin 10 ${strand_field}

    mv ${sample_name}.Aligned.out.bam ${sample_name}.bam

    python ${script_dir}/rnaseq_qc/rnaseq_matrix_plot.py \\
        mapping \\
        ${sample_name}.Log.final.out \\
        ${sample_name}
    """
}


/*
* star bam sort
*/
process star_sortOutput {
    tag "${name}"

    publishDir "${params.outdir}/${params.proj_name}/data/bam/${name}", mode: 'copy',
        saveAs: {filename -> filename.indexOf(".sorted.b*") > 0 ? "${filename}" : null}

    cpus = 8

    input:
    file bam from unsort_bam

    output:
    file "${name}.sorted.bam" into picard_bam, is_bam, bam_assembly, bam_rmats
    file "${name}.sorted.bai" into picard_bam_idx, is_bam_idx

    script:
    name = bam.baseName

    """
    picard ${params.markdup_java_options} FixMateInformation \\
        I=${bam} \\
        O=${name}.sorted.bam \\
        SORT_ORDER=coordinate CREATE_INDEX=true \\
        TMP_DIR=./tmp
    """
}

/*
* picard qc
*/
process picard_qc {
    tag "${name}"

    publishDir "${params.outdir}/${params.proj_name}/result/qc/${name}", mode: 'copy'

    cpus = 8

    input:
    file bam from picard_bam
    file bam_idx from picard_bam_idx
    file ref_flat from ref_flat

    output:
    file "${name}.RnaSeqMetrics.txt" into rnaseq_met, rnaseq_met_qc
    file "${name}.genome_region*" into genome_region_plt, genome_region_plt_qc
    file "${name}.gene_coverage*" into gene_coverage_plt, gene_coverage_plt_qc

    script:
    name = bam.baseName - '.sorted'
    strand_sp = params.picard_lib[params.mrna_strand]
    """
    picard ${params.markdup_java_options} CollectRnaSeqMetrics \\
        I=${bam} \\
        O=${name}.RnaSeqMetrics.txt \\
        REF_FLAT=${ref_flat} \\
        STRAND=${strand_sp}

    python ${script_dir}/rnaseq_qc/rnaseq_matrix_plot.py \\
        rnaseq_matrix \\
        ${name}.RnaSeqMetrics.txt
    """
}


/*
* picard qc
*/

process picard_IS {
    tag "${name}"

    cpus = 8

    input:
    file bam from is_bam
    file bam_idx from is_bam_idx

    output:
    file "${name}.insert_size_metrics.txt" into insert_size, insert_size_qc
    file "${name}.insert_size_histogram.png" into insert_size_plot, insert_size_plot_qc

    script:
    name = bam.baseName - '.sorted'
    """
    picard ${params.markdup_java_options} CollectInsertSizeMetrics \\
        I=${bam} \\
        O=${name}.insert_size_metrics.txt \\
        H=${name}.insert_size_histogram.pdf

    convert -flatten ${name}.insert_size_histogram.pdf ${name}.insert_size_histogram.png
    """
}

/*
* qc report
*/
process qc_report {

    publishDir "${params.outdir}/${params.proj_name}/result/", mode: 'copy',
        saveAs: {filename -> 
            if (filename.indexOf("pping_summary") > 0)  "alignment/${filename}"
            else if (filename.indexOf("gene_coverage") > 0) "qc/${filename}"
            else if (filename.indexOf("genome_region") > 0) "qc/${filename}"
            else if (filename.indexOf("SeqMetrics_summary.csv") > 0) "qc/${filename}"
	        else if (filename == 'report') "../qc_report"
            else null
            }

    input:
    file v_software from software_versions_yaml_qc
    file data_summary from data_summary_qc
    file gc_plot from gc_plot_qc
    file rq_plot from rq_plot_qc
    file rf_plot from rf_plot_qc
    file ('star/*') from star_log_qc.collect()
    file ('star/*') from mapping_plt_qc.collect()
    file ('insert_size/*') from insert_size_qc.collect()
    file ('insert_size/*') from insert_size_plot_qc.collect()
    file ('rnaseq_matrix/*') from rnaseq_met_qc.collect()
    file ('genome_region/*') from genome_region_plt_qc.collect()
    file ('gene_coverage/*') from gene_coverage_plt_qc.collect()

    output:
	file "*.{csv,pdf,png}" into qc_results
    file "report" into qc_report

    when:
    params.qc

    script:
    """
    python ${script_dir}/rnaseq_qc/star_mapping_summary.py \\
        ./star ./

    python ${script_dir}/rnaseq_qc/rnaseq_matrix_summary.py \\
        -d ./rnaseq_matrix -o ./ -m rnaseq

    python ${script_dir}/rnaseq_qc/rnaseq_matrix_summary.py \\
        -d ./insert_size -o ./ -m is

    python ${script_dir}/report/report.py . --report_title "QC report"

    """
    
}

/*
* Transcriptome Assembly
*/
process assembly {

    tag "${name}"

    publishDir "${params.outdir}/${params.proj_name}/result/assembly/detail/${name}", mode: 'copy'

    input:
    file bam from bam_assembly
    file gtf from gtf

    output:
    file "${name}.gtf" into assembly_gtf
    file "cmp2ref.${name}.gtf.tmap" into assembly_tmap

    when:
    params.pipeline

    script:
    name = bam.baseName - '.sorted'
    st_direction = params.stringtie_lib[params.mrna_strand]
    """
    /public/software/stringtie/stringtie-2.0.6.Linux_x86_64/stringtie \\
        -G ${gtf} ${st_direction} \\
        --conservative \\
        -o ${name}.gtf \\
        ${bam} 
        
    gffcompare -r ${gtf} -o cmp2ref ${name}.gtf

    """
}

/*
* Novel Transcripts
*/
process stringtie_merge {

    publishDir "${params.outdir}/${params.proj_name}/result/assembly/", mode: 'copy',
            saveAs: {filename ->
                if (filename == "novel.gtf") filename
                else if (filename == "novel.fa") filename
                else null
            }

    input:
    file "gtf/*" from assembly_gtf.collect()
    file fasta from genome
    file gtf from gtf
    
    output:
    file "novel.gtf" into novel_gtf, novel_anno_gtf
    file "novel.fa" into novel_fasta
    file "analysis.gtf" into kallisto_idx_gtf, quant_gtf, rmats_gtf
    file "assembly_summary.p*" into assembly_unannotated_fig
    file "gene_trans.map" into quant_gene2tr, anno_gene2tr

    when:
    params.pipeline

    script:
    stranded_flag = params.stranded_flag[params.mrna_strand]
    """
    ls gtf/* > gtf.list

    stringtie --merge \\
        -G ${gtf} \\
        -m 200 \\
        -o assembly.gtf \\
        -T 1 \\
        -f 0.05 \\
        gtf.list

    gffcompare \\
        -r ${gtf} \\
        -R assembly.gtf \\
        -o cmp2ref

    python ${script_dir}/assembly/assembly_stats.py stats_combined \\
       --tmap cmp2ref.assembly.gtf.tmap ${stranded_flag} \\
       --outdir .

    python ${script_dir}/assembly/novel_gtf_from_gffcompare.py \\
        --compare-gtf cmp2ref.annotated.gtf ${stranded_flag} \\
        --outfile novel.gtf

    gffread novel.gtf \\
        -g ${fasta} \\
        -w novel.fa    

    cat novel.gtf ${gtf} > analysis.gtf

    python ${script_dir}/assembly/get_gene_to_trans.py \\
        --gff analysis.gtf \\
        --out_dir .
    """
}

/*
* novel gene annotation
*/
process novel_gene_annotation {

    publishDir "${params.outdir}/${params.proj_name}/result/assembly/", mode: 'copy',
            saveAs: {filename ->
                if (filename == "novel.gene.annotation.csv") filename
                else null
            }

    queue 'om'

    cpus = 40

    input:
    file novel_fa from novel_fasta
    file novel_gtf from novel_anno_gtf
    file gene2tr from anno_gene2tr
    file gtf from gtf
    file gene_des from gene_ann
    
    output:
    file "novel.gene.annotation.csv" into novel_gene_out
    file "all.gene.annotation.txt" into quant_gene_ann, diff_gene_ann
    
    script:
    """
    python ${script_dir}/novel_gene/annotate_novel_gene.py \\
        --input-file ${novel_fa} \\
        --gene2tr ${gene2tr}
    
    python ${script_dir}/novel_gene/gene_location.py \\
        --gtf-file ${novel_gtf} \\
        --outfile novel.loc.txt

    python ${script_dir}/novel_gene/gene_location.py \\
        --gtf-file ${gtf} \\
        --outfile known.loc.txt

    python ${script_dir}/utils/merge_files.py \\
        known.loc.txt \\
        ${gene_des} \\
        known.gene.annotation.txt \\
        --na_rep '--' \\
        --method left 


    python ${script_dir}/utils/merge_files.py \\
        novel.loc.txt \\
        novel.annotation.txt \\
        novel.gene.annotation.txt \\
        --na_rep '--' \\
        --method left 

    python ${script_dir}/utils/merge_files.py \\
        known.gene.annotation.txt \\
        novel.gene.annotation.txt \\
        all.gene.annotation.txt \\
        --by_row

    python ${script_dir}/utils/merge_files.py \\
        novel.loc.txt \\
        novel.annotation.detail.txt \\
        novel.gene.annotation.csv \\
        --na_rep '--' \\
        --method left \\
        --out_sep ','
    """
}

/*
* quantification
*/
process mk_kallisto_index {
    tag "KALLISTO index"

    input:
    file gtf from kallisto_idx_gtf
    file fasta from genome
    file fai from genome_fai
    
    output:
    file "transcript.fa" into merged_fa
    file "transcript.fa.kallisto_idx" into kallisto_idx
    
    when:
    params.pipeline

    script:
    """
    gffread ${gtf} \\
        -g ${fasta} \\
        -w transcript.fa

    kallisto index \\
        -i transcript.fa.kallisto_idx \\
        transcript.fa
    """
}


/*
* Kallisto quant
*/
process kallisto {

    tag "${name}"

    publishDir "${params.outdir}/${params.proj_name}/data/kallisto/", mode: 'copy'

    cpus 4
    
    input:
    file reads from quant_reads
    file index from kallisto_idx
    
    output:
    file name into kallisto_out

    when:
    params.pipeline

    script:
    name = reads[0].toString() - '.R1.clean.fq.gz'
    st_direction = params.kallisto_lib[params.mrna_strand]
    """
    kallisto quant ${st_direction} \\
        --threads ${task.cpus} \\
        -i ${index} \\
        --output-dir=${name} \\
        ${reads}   
    """
}

process load_kallisto_results {

    publishDir "${params.outdir}/${params.proj_name}/result/quantification/", mode: 'copy',
        saveAs: {filename -> filename == "deg_input.RData" ? null : "$filename"}

    input:
    file 'kallisto/*' from kallisto_out.collect()
    file sample_group from sample_group
    file gene2tr from quant_gene2tr
    file gene_feature from quant_gene_ann

    output:
    file 'expression_summary' into expression_summary, exp_summary_cor, exp_summary_report
    file 'deg_input.RData' into deg_obj
    file "expression_summary/Sample_correlation_heatmap.png" into cor_plot
    file "expression_summary/PCA_plot.png" into pca_plot

    when:
    params.pipeline
    
    script:
    """   
    /public/software/R/R-3.5.1/executable/bin/Rscript ${script_dir}/quant/mrna_kallisto_to_table.R \\
        --kallisto_dir kallisto \\
        --sample_inf ${sample_group} \\
        --gene2tr ${gene2tr} \\
        --out_dir expression_summary \\
        --gene_feature ${gene_feature}
    """
}

/*
Diff analysis
*/
process rmats_cfg {

    publishDir "${params.outdir}/${params.proj_name}/configure" , mode: 'copy'

    input:
    file sample_group from sample_group
    file contrast_file from contrast_file
    file 'bam/*' from bam_rmats.collect()
    
    output:
    file 'rmats_compare/*' into rmats_compare
    
    when:
    params.pipeline
    
    script:
    """
    python ${script_dir}/utils/make_analysis_compare.py \\
        rmats-sample-files \\
        --bam-dir ./bam \\
        --sample-group ${sample_group} \\
        --out-dir rmats_compare \\
        --contrast ${contrast_file}
    """
}

rmats_compare
    .flatMap()
    .set { diff_compare }

process diff_analysis {
    tag "${compare}"

    publishDir "${params.outdir}/${params.proj_name}/result/quantification/differential_analysis/", mode: 'copy',
        saveAs: {filename -> 
            if (filename.indexOf("png") > 0) null
            else if (filename.indexOf("txt") > 0) null
            else "${filename}"
        }
        
    input:
    file compare from diff_compare
    file deg_obj from deg_obj
    file sample_group from sample_group
    file gf_diff from diff_gene_ann
    
    output:
    file compare into diff_out_go, diff_out_kegg, diff_out_all, diff_tf, diff_out_ppi, diff_out_enrich_plot, pathway_compare, diff_summary, diff_venn
    file "${compare}/${compare}_Volcano_plot.png" into diff_plt
    file "${compare}/${compare}.edgeR.DE_results.txt" into diff_table
    
    when:
    params.pipeline
    
    script:
    """
    mv ${compare} ${compare}_compare

    /public/software/R/R-3.5.1/executable/bin/Rscript ${script_dir}/quant/mrna_diff_edgeR.R \\
        --deg_rdata ${deg_obj} \\
        --compare ${compare} \\
        --sample_inf ${sample_group} \\
        --out_dir ${compare} \\
        --qvalue ${params.exp_diff_pval} \\
        --logfc ${params.exp_lgfc} \\
        --gene_ann ${gf_diff}
    """
}

if (params.chr_size && params.chr_window) {
    /*
    * diff gene distribution
    */
    chr_size = check_ref_exist(params.chr_size, 'chr size file')
    chr_window = check_ref_exist(params.chr_window, 'chr window file')

    process diff_gene_loc {

        tag "${compare}"

        publishDir "${params.outdir}/${params.proj_name}/result/quantification/differential_analysis/${compare}", mode: 'copy'
        
        input:
        file diff_table from diff_table
        file chr_size from chr_size
        file chr_window from chr_window
        
        output:
        file "${compare}_diff_gene_location.*" into all_diff_loc
            
        when:
        params.pipeline

        script:
        compare = diff_table.baseName - '.edgeR.DE_results'
        """
        python ${script_dir}/quant/diff_distribute.py \\
                --diff-file ${diff_table} \\
                --chr-window ${chr_window} \\
                --chr-size ${chr_size} \\
                --pval ${params.exp_diff_pval} \\
                --logfc ${params.exp_lgfc} 
        """
    }

}


/*
* diff summary
*/
process diff_exp_summary {

    publishDir "${params.outdir}/${params.proj_name}/result/", mode: 'copy',
        saveAs: {filename -> 
            if (filename.indexOf("count.txt") > 0) "quantification/expression_summary/${filename}"
            else if (filename.indexOf("tpm.txt") > 0) "quantification/expression_summary/${filename}"
            else if (filename.indexOf("heatmap") > 0) "quantification/expression_summary/${filename}"
            else null
        }

    input:
    file expression_summary from expression_summary
    file "diff/*" from diff_summary.collect()
    file sample_group from sample_group

    output:
    file "Diff_genes_heatmap.{pdf,png}" into diff_heatmap
    file 'Diff_genes*count.txt' into diff_count
    file 'Diff_genes_tpm.txt' into diff_tpm
    file 'cluster_plot' optional true into cluster_plot
    file 'Diff_genes_cluster.txt' optional true into cluster_gene_detail

    when:
    params.pipeline

    script:
    """
    /public/software/R/R-3.5.1/executable/bin/Rscript ${script_dir}/quant/mrna_quant_report.R \\
        --exp_dir ${expression_summary} \\
        --diff_dir diff \\
        --sample_inf ${sample_group} \\
        --out_dir . \\
        --grp_num ${params.min_cluster_grp}
    """
}

/*
* Venn diagram
*/
process diff_gene_venn {

    publishDir "${params.outdir}/${params.proj_name}/result/quantification/expression_summary/", mode: 'copy',
        saveAs: {filename -> 
            if (filename.indexOf(".png") > 0) null
            else filename
        }

    input:
    file "diff/*" from diff_venn.collect()
    file venn_list from venn_list

    output:
    file "venn_plot" into venn_plot_dir
    file "report_plot/*.venn.png" into venn_plot_file
    file "upset_plot" optional true into upset_plot_dir
    file "report_plot/*.upset_plot.png" optional true into upset_plot_file
    
    when:
    params.venn && params.pipeline

    script:
    """
    python ${script_dir}/quant/venn_plot.py \\
        --diff_dir diff \\
        --combination_file ${venn_list} \\
        --out_dir .
    """
}



/*
* enrichment analysis
*/
// go enrichment

def compare2reg = {compare ->
    compare_name = compare.baseName
    return ["${compare_name};ALL", "${compare_name};UP", "${compare_name};DOWN"]
}

diff_out_all
    .flatMap(compare2reg)
    .into { go_compare; kegg_compare; ppi_compare }

process go_analysis {
    tag "${go_compare}"

    errorStrategy 'ignore'

    publishDir "${params.outdir}/${params.proj_name}/result/enrichment/go/", mode: 'copy'

    input:
    file "diff/*" from diff_out_go.collect()
    val go_compare from go_compare
    file go_anno from go_anno
    file gene_length from gene_length
    
    output:
    file "${compare}/${compare}_${reg}_go_enrichment.txt" into go_out
    file "${compare}/DAG/" into dag_plot
    file "${compare}/${compare}_${reg}_go_enrichment_barplot*" into go_barplot

    when:
    params.pipeline    
    
    script:
    (compare, reg) = go_compare.split(';')
    """
    if [ -s diff/${compare}/${compare}.${reg}.edgeR.DE_results.diffgenes.txt ]
    then
        /public/software/R/R-3.5.1/executable/bin/Rscript \\
            ${script_dir}/enrichment/go_analysis.R \\
            --name ${compare}_${reg} \\
            --gene_list diff/${compare}/${compare}.${reg}.edgeR.DE_results.diffgenes.txt \\
            --go_anno ${go_anno} \\
            --gene_length ${gene_length} \\
            --out_dir ${compare}
    fi

    if [ -s ${compare}/${compare}_${reg}_go_enrichment.txt ]
    then
        /public/software/R/R-3.5.1/executable/bin/Rscript \\
            ${script_dir}/enrichment/enrich_bar.R \\
            --enrich_file ${compare}/${compare}_${reg}_go_enrichment.txt
    fi
    """
}


//KEGG enrichment
process kegg_analysis {

    tag "${compare_reg}"

    errorStrategy 'ignore'

    conda "/public/software/miniconda3/envs/kobas/"

    queue 'om'
    
    publishDir "${params.outdir}/${params.proj_name}/result/enrichment/kegg/", mode: 'copy',
        saveAs: {filename -> filename.indexOf("kegg_enrichment") > 0 ? "$filename" : null}   

    input:
    val compare_reg from kegg_compare
    file "diff/*" from diff_out_kegg.collect()
    file kegg_anno from kegg_anno
    
    output:
    file "${compare}/${compare}_${reg}_kegg_enrichment.txt" into kegg_out, kegg_pathway_input
    file "${compare}.${reg}.blasttab" into kegg_pathway_blast
    file "${compare}/${compare}_${reg}_kegg_enrichment_barplot*" into kegg_barplot

    when:
    params.pipeline && params.kegg_abbr
    
    script:
    (compare, reg) = compare_reg.split(';')
    """
    if [ -s diff/${compare}/${compare}.${reg}.edgeR.DE_results.diffgenes.txt ]
    then
        python ${script_dir}/utils/extract_info_by_id.py \\
            --id diff/${compare}/${compare}.${reg}.edgeR.DE_results.diffgenes.txt \\
            --table ${kegg_anno} \\
            --output ${compare}.${reg}.blasttab
    fi

    if [ -s ${compare}.${reg}.blasttab ]
    then
        mkdir -p ${compare}
        kobas-run \\
            -i ${compare}.${reg}.blasttab \\
            -t blastout:tab \\
            -s ${params.kegg_abbr} \\
            -b ${kegg_backgroud} \\
            -d K \\
            -o ${compare}/${compare}_${reg}_kegg_enrichment.txt
    fi

    if [ -s ${compare}/${compare}_${reg}_kegg_enrichment.txt ]
    then
        python ${script_dir}/enrichment/treat_kegg_table.py \\
            --kegg ${compare}/${compare}_${reg}_kegg_enrichment.txt

        /public/software/R/R-3.5.1/executable/bin/Rscript ${script_dir}/enrichment/enrich_bar.R \\
            --enrich_file ${compare}/${compare}_${reg}_kegg_enrichment.txt
    fi
    """
}


/*
* analysis report
*/

process pipe_report {

    publishDir "${params.outdir}/${params.proj_name}/result/", mode: 'copy',
        saveAs: {filename -> 
            if (filename.indexOf("pping_summary") > 0)  "alignment/${filename}"
            else if (filename.indexOf("gene_coverage") > 0) "qc/${filename}"
            else if (filename.indexOf("genome_region") > 0) "qc/${filename}"
            else if (filename.indexOf("SeqMetrics_summary.csv") > 0) "qc/${filename}"
	        else if (filename == 'report') "../analysis_report"
            else null
            }

    input:
    file v_software from software_versions_yaml
    file data_summary from data_summary
    file gc_plot from gc_plot
    file rq_plot from rq_plot
    file rf_plot from rf_plot
    file ('star/*') from star_log.collect()
    file ('star/*') from mapping_plt.collect()
    file ('rnaseq_matrix/*') from rnaseq_met.collect()
    file ('genome_region/*') from genome_region_plt.collect()
    file ('gene_coverage/*') from gene_coverage_plt.collect()
    file assembly_unannotated_fig
    file ('assembly_tmap/*') from assembly_tmap.collect()
    file expression_summary from exp_summary_report
    file pca_plot from pca_plot
    file cor_plot from cor_plot
    file ('volcano/*') from diff_plt.collect()
    file ('diff_loc/*') from all_diff_loc.collect()
    file ('heatmap/*') from diff_heatmap
    file cluster_plot from cluster_plot.ifEmpty { null }
    file ('venn_plot/*') from venn_plot_file
    file ('upset_plot/*') from upset_plot_file.ifEmpty { null }
    file ('go_barplot/*') from go_barplot.collect()
    file ('kegg_barplot/*') from kegg_barplot.collect()
    


    output:
	file "*.{csv,pdf,png}" into analysis_results
    file "report" into analysis_report

    when:
    params.pipeline

    script:
    stranded_flag = params.stranded_flag[params.mrna_strand]
    """
    python ${script_dir}/rnaseq_qc/star_mapping_summary.py \\
        ./star ./

    python ${script_dir}/rnaseq_qc/rnaseq_matrix_summary.py \\
        -d ./rnaseq_matrix -o ./ -m rnaseq

    python ${script_dir}/assembly/assembly_stats.py stats_all \\
        --tmap_dir assembly_tmap ${stranded_flag} \\
        --out_dir .

    python ${script_dir}/report/report.py . --report_title "转录组分析报告"

    """
    
}
