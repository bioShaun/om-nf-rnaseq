#!/usr/bin/env nextflow

// TODO LIST
// - different libtype choice
// - optional analysis flag: rmats / snp / ...

def helpMessage() {
    log.info """

    A RNAseq analysis pipeline

    Usage:

    ==================================================================

    References
      --hisat2_index                Path to reference hisat index
      --genome                      Path to reference genome fa
      --gtf                         Path to reference gtf file
      --bed12                       Path to reference bed file
    
    Mandatory arguments:
      --proj_name                   Project name for output name
      --reads                       Path to input data (must be surrounded with quotes)
      --sample_group                sample vs group file, tab seperated
      --outdir                      Path to analysis results 

    Optional:
      --skip_rseqc                  Skip rseqc anaysis
      --pipeline                    run whole rnaseq pipeline

    """.stripIndent()
}

// workflow internal path&files
script_dir = file("${baseDir}/script/")
cfg_dir = file("${baseDir}/cfg/")
ch_mdsplot_header = Channel.fromPath("${baseDir}/assets/mdsplot_header.txt")
ch_heatmap_header = Channel.fromPath("${baseDir}/assets/heatmap_header.txt")

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
/*
 * SET UP CONFIGURATION VARIABLES
*/



// parameters
params.proj_name = 'test_project'
params.pathway_db = '/public/database/kegg_html/'
params.outdir = false
params.reads = false
params.hisat2_index = false
params.gtf = false
params.bed12 = false
params.genome = false
params.skip_rseqc = false
params.skip_dupradar = false
params.forward_stranded = false
params.reverse_stranded = false
params.unstranded = true
params.pipeline = false
params.contrast = false
params.exp_diff_pval = 0.05
params.exp_lgfc = 1
params.go_anno = false
params.skip_kegg = false
params.kegg_anno = false
params.kegg_pathway = false
params.gene_length = false
params.kegg_abbr = false
params.kegg_backgroud = false

forward_stranded = params.forward_stranded
reverse_stranded = params.reverse_stranded
unstranded = params.unstranded
params.sample_group = false

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
// hisat index
genome = check_ref_exist(params.genome, 'genome fasta')
gtf = check_ref_exist(params.gtf, 'gtf file')
bed12 = check_ref_exist(params.bed12, 'bed file')
hisat2_index = Channel
                    .fromPath("${params.hisat2_index}/*.ht2*")
                    .ifEmpty { exit 1, "HISAT2 index not found: ${params.hisat2_index}" }
alignment_splicesites = file("${params.hisat2_index}/${genome.baseName}.hisat2_splice_sites.txt")



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
        
        script:
        """
        python ${script_dir}/utils/make_analysis_compare.py \\
            generate_contrast \\
            --sample-group ${sample_group} \\
            --contrast-file contrast.ini
        """
    }
}



// software version
process get_software_versions {

    output:
    file 'software_versions_mqc.yaml' into software_versions_yaml

    script:
    """
    echo $workflow.nextflow.version &> v_nextflow.txt
    fastp --version &> v_fastp.txt
    hisat2 --version &> v_hisat2.txt
    samtools --version &> v_samtools.txt
    stringtie --version &> v_stringtie.txt
    featureCounts -v &> v_featurecounts.txt
    picard MarkDuplicates --version &> v_markduplicates.txt  || true
    multiqc --version &> v_multiqc.txt
    python ${script_dir}/utils/scrape_software_versions.py &> software_versions_mqc.yaml
    """
}

/*
 * Fastp
 */

process fastp {

    tag "${name}"

    publishDir "${params.outdir}/${params.proj_name}/data/fastp_trimmed_reads/${name}", mode: 'copy'

    input:
    set name, file(reads) from raw_fq_files

    output:
    file "*trimmed.R*.fq.gz" into trimmed_reads, kallisto_reads
    file "${name}.fastp.json" into fastp_json_analysis, fastp_json_report
    file "${name}.fastp.html" into fastp_html

    cpus = 4

    script:
    """
    fastp \\
        --in1 ${reads[0]} \\
        --in2 ${reads[1]} \\
        --out1 ${name}.trimmed.R1.fq.gz \\
        --out2 ${name}.trimmed.R2.fq.gz \\
        --json ${name}.fastp.json \\
        --html ${name}.fastp.html    
    """
}  


/*
* fq quality control
*/
process reads_qc_summary {

    publishDir "${params.outdir}/${params.proj_name}/result/reads_qc/", mode: 'copy'

    input:
    file 'fastp_json/*' from fastp_json_analysis.collect()
    
    output:
    file "data.summary.csv"
    file "reads_gc"
    file "reads_quality"
    file "reads_filter"
    
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

process hisat_mapping {
    tag "${sample_name}"

    publishDir "${params.outdir}/${params.proj_name}/result/alignment/${sample_name}", mode: 'copy',
        saveAs: {filename -> filename.indexOf("hisat2_summary.txt") > 0 ? "$filename" : null}

    input:
    file reads from trimmed_reads
    file index from hisat2_index.collect()
    
    output:
    file "${sample_name}.bam" into unsort_bam
    file "${sample_name}.hisat2_summary.txt" into alignment_logs

    cpus = 8
    
    script:
    sample_name = reads[0].toString() - '.trimmed.R1.fq.gz'
    index_base = index[0].toString() - ~/.\d.ht2l?/
    """
    hisat2 -x ${index_base} \\
            -1 ${reads[0]} \\
            -2 ${reads[1]} \\
            --known-splicesite-infile ${alignment_splicesites} \\
            --no-mixed \\
            --no-discordant \\
            -p ${task.cpus} \\
            --met-stderr \\
            --new-summary \\
            --summary-file ${sample_name}.hisat2_summary.txt \\
            --rg-id ${sample_name} \\
            --rg SM:${sample_name}\\
            --rg LB:${sample_name} \\
            --rg PI:350 \\
            --rg PL:Illumina \\
            --rg CN:TCuni \\
            | samtools view -bS -F 4 -F 8 -F 256 - > ${sample_name}.bam           
    """
}


/*
* Hisat2 bam sort
*/
process hisat2_sortOutput {
    tag "${name}"

    publishDir "${params.outdir}/${params.proj_name}/data/bam/${name}", mode: 'copy',
        saveAs: {filename -> filename.indexOf(".sorted.bam") > 0 ? "$filename" : null}

    input:
    file hisat2_bam from unsort_bam

    output:
    file "${name}.sorted.bam" into bam_count, bam_scallop, bam_for_genebody, bam_saturation, bam_rmats, bam_snp, bam_markduplicates, bam_stringtieFPKM, bam_featurecounts
    file "${name}.sorted.bam.bai" into bam_for_genebody_idx
    file "${name}.fixmate.bam"

    cpus = 8

    script:
    name = hisat2_bam.baseName
    """
    samtools fixmate --threads ${task.cpus} \\
        -m ${hisat2_bam} ${name}.fixmate.bam

    samtools sort \\
        ${name}.fixmate.bam \\
        --threads ${task.cpus} \\
        -o ${name}.sorted.bam
    
    samtools index ${name}.sorted.bam
    """
}


/*
* Rseqc create BigWig coverage
*/
process rseqc {
    tag "${name}"

    publishDir "${params.outdir}/${params.proj_name}/result/rnaseq_qc/" , mode: 'copy',
        saveAs: {filename ->
                 if (filename.indexOf("infer_experiment.txt") > 0)              "infer_experiment/$filename"
            //else if (filename.indexOf("eRPKM.xls") > 0)         "RPKM_saturation/rpkm/$filename"
            //else if (filename.indexOf("rawCount.xls") > 0)      "RPKM_saturation/counts/$filename"
            //else if (filename.indexOf("saturation.pdf") > 0)    "RPKM_saturation/$filename"
            //else if (filename.indexOf("saturation.r") > 0)      "RPKM_saturation/rscripts/$filename"
            //else if (filename.indexOf("junction_plot.r") > 0)                   "junction_annotation/rscripts/$filename"
            else if (filename.indexOf("junction.xls") > 0)                      "junction_annotation/data/$filename"
            else if (filename.indexOf("splice_events.pdf") > 0)                 "junction_annotation/events/$filename"
            else if (filename.indexOf("splice_junction.pdf") > 0)               "junction_annotation/junctions/$filename"
            else if (filename.indexOf("junctionSaturation_plot.pdf") > 0)       "junction_saturation/$filename"
            //else if (filename.indexOf("junctionSaturation_plot.r") > 0)         "junction_saturation/rscripts/$filename"
            else if (filename.indexOf("inner_distance.txt") > 0)                "inner_distance/$filename"
            else if (filename.indexOf("inner_distance_freq.txt") > 0)           "inner_distance/data/$filename"
            //else if (filename.indexOf("inner_distance_plot.r") > 0)             "inner_distance/rscripts/$filename"
            else if (filename.indexOf("inner_distance_plot.pdf") > 0)           "inner_distance/plots/$filename"
            else null
        }

    when:
    !params.skip_qc && !params.skip_rseqc

    input:
    file bam from bam_for_genebody
    file bam_idx from bam_for_genebody_idx
    file bed12 from bed12

    output:
    file "*.bigwig" into bigwig_for_genebody
    file "*.{txt,pdf,r,xls}" into rseqc_results

    cpus 8

    script:
    name = bam.baseName - '.sorted'
    """    
    bamCoverage -b ${bam} -p ${task.cpus} -o ${name}.bigwig
    infer_experiment.py -i ${bam} -r ${bed12} > ${name}.infer_experiment.txt
    junction_annotation.py -i ${bam} -o ${name}.rseqc -r ${bed12}
    junction_saturation.py -i ${bam} -o ${name}.rseqc -r ${bed12} 2> ${name}.junction_annotation_log.txt
    RPKM_saturation.py -i ${bam} -o ${name} -r ${bed12}   
    inner_distance.py -i ${bam} -o ${name}.rseqc -r ${bed12}
    """
}


/*
 * Rseqc genebody coverage
 */
process genebody_coverage {
    tag "${name}"

    publishDir "${params.outdir}/${params.proj_name}/result/rnaseq_qc/genebody_cov" , mode: 'copy'

    when:
    !params.skip_qc && !params.skip_rseqc

    input:
    file bigwig from bigwig_for_genebody
    file bed12 from bed12

    output:
    file "${name}.geneBodyCoverage.txt" into genebody_coverage_results

    cpus 8

    script:
    name = bigwig.baseName
    """
    geneBody_coverage2.py \\
        -i ${bigwig} \\
        -o ${name} \\
        -r ${bed12}

    python ${script_dir}/rnaseq_qc/format_genecov.py ${name}.geneBodyCoverage.txt
    """
}


/*
 * Rseqc genebody coverage plot
 */
/* process genebody_coverage_plot {

    publishDir "${params.outdir}/${params.proj_name}/result/rnaseq_qc" , mode: 'copy'

    input:
    file "genebody_coverage/*" from genebody_coverage_results.collect()

    output:
    file "genebody_coverage" into genebody_coverage


    script:
    """
    /public/software/R/R-3.5.1/executable/bin/Rscript ${script_dir}/rnaseq_qc/genebody_cov.R \\
        --cov_dir genebody_coverage
    """
} */


/*
 * Rseqc gene exp saturation analysis
 */
/* process saturation_plot {

    publishDir "${params.outdir}/${params.proj_name}/result/rnaseq_qc/" , mode: 'copy'

    input:
    file 'sequencing_saturation/*' from count_saturation.collect()

    output:
    file "sequencing_saturation" into sequencing_saturation

    script:
    """
    /public/software/R/R-3.5.1/executable/bin/Rscript ${script_dir}/rnaseq_qc/seq_saturation.R \\
        --saturation_count sequencing_saturation
    """
} */

process markDuplicates {
    tag "${name}"
    publishDir "${params.outdir}/${params.proj_name}/data/bam/${name}", mode: 'copy'

    when:
    !params.skip_qc && !params.skip_dupradar

    input:
    file bam from bam_markduplicates

    output:
    file "${name}.markDups.bam" into bam_md
    file "${name}.markDups_metrics.txt" into picard_results
    file "${name}.markDups.bam.bai"

    script:
    name = bam.baseName - '.sorted'
    markdup_java_options = (task.memory.toGiga() > 8) ? ${params.markdup_java_options} : "\"-Xms" +  (task.memory.toGiga() / 2 )+"g "+ "-Xmx" + (task.memory.toGiga() - 1)+ "g\""

    """
    picard ${markdup_java_options} MarkDuplicates \\
        INPUT=${bam} \\
        OUTPUT=${name}.markDups.bam \\
        METRICS_FILE=${name}.markDups_metrics.txt \\
        REMOVE_DUPLICATES=false \\
        ASSUME_SORTED=true \\
        PROGRAM_RECORD_ID='null' \\
        VALIDATION_STRINGENCY=LENIENT
    samtools index ${name}.markDups.bam
    """
}


process dupradar {
    tag "${name}"
    publishDir "${params.outdir}/${params.proj_name}/result/rnaseq_qc/dupradar", mode: 'copy',
        saveAs: {filename ->
            if (filename.indexOf("_duprateExpDens.pdf") > 0) "scatter_plots/$filename"
            else if (filename.indexOf("_duprateExpBoxplot.pdf") > 0) "box_plots/$filename"
            else if (filename.indexOf("_expressionHist.pdf") > 0) "histograms/$filename"
            else if (filename.indexOf("_dupMatrix.txt") > 0) "gene_data/$filename"
            else if (filename.indexOf("_duprateExpDensCurve.txt") > 0) "scatter_curve_data/$filename"
            else if (filename.indexOf("_intercept_slope.txt") > 0) "intercepts_slopes/$filename"
            else "$filename"
        }

    when:
    !params.skip_qc && !params.skip_dupradar

    input:
    file bam_md
    file gtf from gtf

    output:
    file "*.{pdf,txt}" into dupradar_results

    script: 
    name = bam_md.baseName - '.markDups'
    def dupradar_direction = 0
    if (forward_stranded && !unstranded) {
        dupradar_direction = 1
    } else if (reverse_stranded && !unstranded){
        dupradar_direction = 2
    }
    """
    /public/software/R/R-3.5.1/executable/bin/Rscript ${script_dir}/rnaseq_qc/dupRadar.r ${bam_md} ${gtf} ${dupradar_direction} paired ${task.cpus}
    """
}


process featureCounts {
    label 'low_memory'
    tag "${name}"
    publishDir "${params.outdir}/${params.proj_name}/result/rnaseq_qc/featureCounts", mode: 'copy',
        saveAs: {filename ->
            if (filename.indexOf("biotype_counts") > 0) "biotype_counts/$filename"
            else if (filename.indexOf("_gene.featureCounts.txt.summary") > 0) "gene_count_summaries/$filename"
            else if (filename.indexOf("_gene.featureCounts.txt") > 0) "gene_counts/$filename"
            else "$filename"
        }

    when:
    !params.skip_qc && !params.skip_expCorr  

    input:
    file bam_featurecounts
    file gtf from gtf

    output:
    file "${name}_gene.featureCounts.txt" into geneCounts, featureCounts_to_merge
    file "${name}_gene.featureCounts.txt.summary" into featureCounts_logs
    
    script:
    def featureCounts_direction = 0
    def extraAttributes = params.fcExtraAttributes ? "--extraAttributes ${params.fcExtraAttributes}" : ''
    if (forward_stranded && !unstranded) {
        featureCounts_direction = 1
    } else if (reverse_stranded && !unstranded){
        featureCounts_direction = 2
    }
    // Try to get real sample name
    name = bam_featurecounts.baseName - '.sorted'
    """
    featureCounts -a ${gtf} -g ${params.fcGroupFeatures} -o ${name}_gene.featureCounts.txt ${extraAttributes} -p -s ${featureCounts_direction} ${bam_featurecounts}
    """
}



process stringtieFPKM {
    tag "${name}"
    publishDir "${params.outdir}/${params.proj_name}/result/stringtie/${name}", mode: 'copy',
        saveAs: {filename ->
            if (filename.indexOf(".gtf") > 0) "${filename}"
            else null
        }

    input:
    file bam_stringtieFPKM
    file gtf from gtf

    output:
    file "${name}.gtf" into assembly_gtf
    file "${name}.gffcompare.${name}.gtf.tmap" into stringtie_log

    script:
    name = bam_stringtieFPKM.baseName - '.sorted'
    def st_direction = ''
    if (forward_stranded && !unstranded){
        st_direction = "--library_type second"
    } else if (reverse_stranded && !unstranded){
        st_direction = "--library_type first"
    }
    """
    scallop \\
        -i ${bam_stringtieFPKM} \\
        ${st_direction} \\
        -o ${name}.gtf

    gffcompare -r ${gtf} -o ${name}.gffcompare ${name}.gtf

    """
}


process sample_correlation {

    publishDir "${params.outdir}/${params.proj_name}/result/rnaseq_qc/sample_correlation", mode: 'copy'

    when:
    !params.skip_qc && !params.skip_expCorr

    input:
    file input_files from geneCounts.collect()
    val num_bams from bam_count.count()
    file mdsplot_header from ch_mdsplot_header
    file heatmap_header from ch_heatmap_header

    output:
    file "*.{txt,pdf,csv}" into sample_correlation_results

    when:
    num_bams > 2 && (!params.sampleLevel)

    script: // This script is bundled with the pipeline, in nfcore/rnaseq/bin/
    """
    /public/software/R/R-3.5.1/executable/bin/Rscript ${script_dir}/rnaseq_qc/edgeR_heatmap_MDS.r ${input_files}
    """
}



process stringtie_merge {

    publishDir "${params.outdir}/${params.proj_name}", mode: 'copy',
            saveAs: {filename ->
                if (filename == "novel.gtf") "result/stringtie/$filename"
                else if (filename == "assembly.gtf") "result/stringtie/$filename"
                else if (filename == "transcripts.gtf") "data/ref/$filename"
                else null
            }

    when:
    params.pipeline

    input:
    file "gtf/*" from assembly_gtf.collect()
    file gtf from gtf
    
    output:
    file "transcripts.gtf" into kallisto_idx_gtf, quant_gtf, rmats_gtf
    file "novel.gtf" into novel_gtf, novel_anno_gtf
    file "assembly.gtf"
    file "gene_trans.map" into gene2tr_quant, gene2tr_anno

    script:
    """
    ls gtf/* > gtf.list

    stringtie --merge \\
        -G ${gtf} \\
        -m 200 \\
        -o assembly.gtf \\
        -T 1 \\
        -f 0.1 \\
        gtf.list

    gffcompare \\
        -r ${gtf} \\
        -R assembly.gtf \\
        -o cmp2ref

    python ${script_dir}/assembly/novel_gtf_from_gffcompare.py \\
        --compare-gtf cmp2ref.annotated.gtf \\
        --outfile novel.gtf

    cat ${gtf} novel.gtf > transcripts.gtf

    python ${script_dir}/assembly/get_gene_to_trans.py \\
        --gff transcripts.gtf \\
        --out_dir .
    """
}


process mk_kallisto_index {
    tag "KALLISTO index"

    publishDir "${params.outdir}/${params.proj_name}", mode: 'copy',
        saveAs: {filename ->
            if (filename == "novel.fa") "result/stringtie/$filename"
            else if (filename == "transcripts.fa.kallisto_idx") "data/ref/$filename"
            else if (filename == "transcripts.fa") "data/ref/$filename"
            else null
        }

    when:
    params.pipeline

    input:
    file gtf from kallisto_idx_gtf
    file novel_gtf from novel_gtf
    file fasta from genome
    
    output:
    file "transcripts.fa" into merged_fa
    file "novel.fa" into novel_fa
    file "transcripts.fa.kallisto_idx" into kallisto_idx
    
    script:
    """
    gffread ${gtf} \\
        -g ${fasta} \\
        -w transcripts.fa

    gffread ${novel_gtf} \\
        -g ${fasta} \\
        -w novel.fa    

    kallisto index \\
        -i transcripts.fa.kallisto_idx \\
        transcripts.fa
    """
}

process kallisto {

    tag "${name}"

    publishDir "${params.outdir}/${params.proj_name}/data/kallisto/", mode: 'copy'

    when:
    params.pipeline    

    input:
    file reads from kallisto_reads
    file index from kallisto_idx
    
    output:
    file name into kallisto_out

    cpus 4
    
    script:
    name = reads[0].toString() - '.trimmed.R1.fq.gz'
    def st_direction = ''
    if (forward_stranded && !unstranded){
        st_direction = "--fr-stranded"
    } else if (reverse_stranded && !unstranded){
        st_direction = "--rf-stranded"
    }
    """
    kallisto quant \
        --threads ${task.cpus} \
        -i ${index} \
        --output-dir=${name} \
        ${st_direction} \
        ${reads}   
    """
}


process rmats_cfg {

    publishDir "${params.outdir}/${params.proj_name}/configure" , mode: 'copy'

    input:
    file sample_group from sample_group
    file contrast_file from contrast_file
    file 'bam/*' from bam_rmats.collect()
    
    output:
    file 'rmats_compare/*' into rmats_compare
    
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
    .into { flat_rmats_compare; diff_compare }


process load_kallisto_results {

    publishDir "${params.outdir}/${params.proj_name}/result/quantification/", mode: 'copy',
        saveAs: {filename -> filename == "deg_input.RData" ? null : "$filename"}

    when:
    params.pipeline    

    input:
    file 'kallisto/*' from kallisto_out.collect()
    file sample_group from sample_group
    file gene2tr from gene2tr_anno
    
    output:
    file 'expression_summary' into expression_summary
    file 'deg_input.RData' into deg_obj
    
    script:
    """   
    /public/software/R/R-3.5.1/executable/bin/Rscript ${script_dir}/quant/kallisto_to_table.R \\
        --kallisto_dir kallisto \\
        --sample_inf ${sample_group} \\
        --gene2tr ${gene2tr} \\
        --out_dir expression_summary
    """
}


/*
* Differential analysis
*/
process diff_analysis {
    tag "${compare}"

    publishDir "${params.outdir}/${params.proj_name}/result/quantification/differential_analysis/", mode: 'copy'

    when:
    params.pipeline 

    input:
    file compare from diff_compare
    file deg_obj from deg_obj
    file sample_group from sample_group
    
    output:
    file compare into diff_out_go, diff_out_kegg, diff_out_all, diff_tf, diff_out_ppi, diff_out_enrich_plot, pathway_compare, diff_summary
    
    script:
    """
    mv ${compare} ${compare}_compare

    /public/software/R/R-3.5.1/executable/bin/Rscript ${script_dir}/quant/diff_edgeR.R \\
        --deg_rdata ${deg_obj} \\
        --compare ${compare} \\
        --sample_inf ${sample_group} \\
        --out_dir ${compare} \\
        --qvalue ${params.exp_diff_pval} \\
        --logfc ${params.exp_lgfc}
    """
}


/*
* diff summary
*/
process diff_exp_summary {

    publishDir "${params.outdir}/${params.proj_name}/result/quantification/expression_summary", mode: 'copy'

    when:
    params.pipeline 

    input:
    file expression_summary from expression_summary
    file "diff/*" from diff_summary.collect()
    file sample_group from sample_group

    output:
    file "Diff.gene*" into diff_summary_out
    file "cluster_data"
    file "diff.genes.*.csv" into diff_genes_stats
    
    script:
    if (params.contrast){
        contrast_flag = "True"
    } else {
        contrast_flag = "True"
    }

    """
    /public/software/R/R-3.5.1/executable/bin/Rscript ${script_dir}/quant/quant_report.R \\
        --exp_dir ${expression_summary} \\
        --diff_dir diff \\
        --sample_inf ${sample_group} \\
        --out_dir .

    python ${script_dir}/quant/diff_stats.py \\
        --diff-dir diff \\
        --group-inf ${sample_group} \\
        --contrasts ${contrast_flag}  \\
        --outdir .
    """
}



/*
* enrichment analysis
*/
// go enrichment

def compare2reg = {compare ->
    compare_name = compare.baseName
    (up_name, down_name) = "${compare.baseName}".split('_vs_')
    return ["${compare_name};ALL", "${compare_name};${up_name}-UP", "${compare_name};${down_name}-UP"]
}

diff_out_all
    .flatMap(compare2reg)
    .into { go_compare; kegg_compare; ppi_compare }


process go_analysis {
    tag "${go_compare}"

    publishDir "${params.outdir}/${params.proj_name}/result/enrichment/go/", mode: 'copy'

    when:
    params.pipeline 

    input:
    file "diff/*" from diff_out_go.collect()
    val go_compare from go_compare
    file go_anno from file(params.go_anno)
    file gene_length from file(params.gene_length)
    
    output:
    file "${compare}/${compare}.${reg}.go.enrichment.txt" into go_out
    file "${compare}/DAG/*" into dag_plot
    file "${compare}/${compare}.${reg}.go.enrichment.barplot*" into go_barplot
    
    script:
    (compare, reg) = go_compare.split(';')
    """
    /public/software/R/R-3.5.1/executable/bin/Rscript ${script_dir}/enrichment/go_analysis.R \\
        --name ${compare}.${reg} \\
        --gene_list diff/${compare}/${compare}.${reg}.edgeR.DE_results.diffgenes.txt \\
        --go_anno ${go_anno} \\
        --gene_length ${gene_length} \\
        --out_dir ${compare}

    /public/software/R/R-3.5.1/executable/bin/Rscript ${script_dir}/enrichment/enrich_bar.R \\
        --enrich_file ${compare}/${compare}.${reg}.go.enrichment.txt
    """
}


process go_report {

    when:
    params.pipeline     

    input:
    file "go/*" from go_out.collect()

    output:
    file "top.go.enrich.json" into go_report_data

    script:
    """
    python ${script_dir}/enrichment/enrich_multiqc_data.py ./go go .
    """

}


//KEGG enrichment
process kegg_analysis {

    tag "${compare_reg}"

    conda '/public/software/miniconda3/envs/kobas/'

    publishDir "${params.outdir}/${params.proj_name}/result/enrichment/kegg/", mode: 'copy',
        saveAs: {filename -> filename.indexOf("kegg.enrichment") > 0 ? "$filename" : null}

    when:
    params.pipeline && !params.skip_kegg

    input:
    val compare_reg from kegg_compare
    file "diff/*" from diff_out_kegg.collect()
    file kegg_anno from file(params.kegg_anno)
    
    output:
    file "${compare}/${compare}.${reg}.kegg.enrichment.txt" into kegg_out, kegg_pathway_input, kegg_out_for_report
    file "${compare}.${reg}.blasttab" into kegg_pathway_blast
    file "${compare}/${compare}.${reg}.kegg.enrichment.barplot*"
    
    script:
    (compare, reg) = compare_reg.split(';')
    """
    python ${script_dir}/utils/extract_info_by_id.py \\
        --id diff/${compare}/${compare}.${reg}.edgeR.DE_results.diffgenes.txt \\
        --table ${kegg_anno} \\
        --output ${compare}.${reg}.blasttab 

    mkdir -p ${compare}
    
    kobas-run \\
        -i ${compare}.${reg}.blasttab \\
        -t blastout:tab \\
        -s ${params.kegg_abbr} \\
        -b ${kegg_backgroud} \\
        -d K \\
        -o ${compare}/${compare}.${reg}.kegg.enrichment.txt

    python ${script_dir}/enrichment/treat_kegg_table.py \\
        --kegg ${compare}/${compare}.${reg}.kegg.enrichment.txt

    /public/software/R/R-3.5.1/executable/bin/Rscript ${script_dir}/enrichment/enrich_bar.R \\
        --enrich_file ${compare}/${compare}.${reg}.kegg.enrichment.txt
    """
}


process kegg_report {

    when:
    params.pipeline     

    input:
    file "kegg/*" from kegg_out_for_report.collect()

    output:
    file "top.kegg.enrich.json" into kegg_report_data

    script:
    """
    python ${script_dir}/enrichment/enrich_multiqc_data.py ./kegg kegg .
    """

}


/*
* KEGG pathway
*/
if (params.kegg_pathway) {
    process kegg_pathway {
        tag "${compare_reg}"

        publishDir "${params.outdir}/${params.proj_name}/result/enrichment/kegg/${compare}", mode: 'copy'

        when:
        params.pipeline 

        input:
        file kegg_out from kegg_pathway_input
        file blast_out from file(params.kegg_anno)
        file "diff/*" from pathway_compare.collect()
        
        output:
        file "${compare_reg}.pathway" into kegg_pathway_out
        
        script:
        compare_reg = kegg_out.baseName - '.kegg.enrichment'
        compare = kegg_out.baseName - ~/.(\w+-UP|ALL).kegg.enrichment$/
        """
        mkdir ${compare_reg}.pathway
        python ${script_dir}/enrichment/pathway_html.py \\
            --pathway-db ${params.pathway_db} \\
            --blasttab ${blast_out} \\
            --kegg-table ${kegg_out} \\
            --diff-table diff/${compare}/${compare}.edgeR.DE_results.txt \\
            --outdir ${compare_reg}.pathway \\
            --kegg-abbr ${params.kegg_abbr} \\
        """
    }
}



/*
* Report
*/
process multiqc {

    publishDir "${params.outdir}/${params.proj_name}/report/", mode: 'copy'

    input:
    file cfg from file("${cfg_dir}/multiqc.cfg.yml")
    file ('fastp/*') from fastp_json_report.collect().ifEmpty([])
    file ('alignment/*') from alignment_logs.collect().ifEmpty([])
    file ('rseqc/*') from rseqc_results.collect().ifEmpty([])
    file ('rseqc/*') from genebody_coverage_results.collect().ifEmpty([])
    file ('dupradar/*') from dupradar_results.collect().ifEmpty([])
    file ('stringtie/*') from stringtie_log.collect().ifEmpty([])
    file ('diff/*') from diff_genes_stats.collect().ifEmpty([])
    file ('go/*') from go_report_data.collect().ifEmpty([])
    file ('kegg/*') from kegg_report_data.collect().ifEmpty([])
    file ('sample_correlation_results/*') from sample_correlation_results.collect().ifEmpty([]) // If the Edge-R is not run create an Empty array
    file ('software_versions/*') from software_versions_yaml
    
    output:
    file "rnaseq_report.html" into multiqc_report
    file "rnaseq_report_data"
    
    script:
    """
    multiqc -t oms -c ${cfg} -n rnaseq_report .
    """
}


