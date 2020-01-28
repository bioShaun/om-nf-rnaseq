// TODO LIST
// - different libtype choice
// - optional analysis flag: rmats / snp / ...

def helpMessage() {
    log.info """

    A RNAseq analysis pipeline

    Usage:

    ==================================================================

    References
      --genome                      Path to reference genome fa
      --gtf                         Path to reference gtf file
      --split_gtf_dir               Path to gtf directory split by gene_type
      --bed12                       Path to reference bed file
      --ref_flat                    Path to reference ref flat
      --star_index                  Path to reference star index
      --rrna_index                  Path to rRNA bowtie2 index
      --gene_ann                    Path to gene annotation
    
    Mandatory arguments:
      --proj_name                   Project name for output name
      --reads                       Path to input data (must be surrounded with quotes)
      --sample_group                Sample vs Group file, tab seperated
      --outdir                      Path to analysis results 

    Optional:
      --rm_rrna                     Remove rRNA sequence in reads
      --qc                          Run qc analysis
      --pipeline                    Run whole rnaseq pipeline

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

/*
 * SET UP CONFIGURATION VARIABLES
*/

// parameters
params.proj_name = 'test_project'
params.pathway_db = '/public/database/kegg_html/'
params.outdir = false
params.reads = false
params.star_index = false
params.gtf = false
params.split_gtf_dir = false
params.ref_flat = false
params.bed12 = false
params.genome = false
params.strand = 'unstranded'
params.contrast = false
params.diff_pval = 0.05
params.diff_lgfc = 1
params.go_anno = false
params.skip_kegg = false
params.kegg_anno = false
params.kegg_pathway = false
params.gene_length = false
params.kegg_abbr = false
params.kegg_backgroud = false
params.sample_group = false
params.rm_rrna = false
params.rrna_index = false
params.gene_ann = false

// pipeline control parameters
params.pipeline = false
params.qc = false
params.quant = false

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

// check ref
genome = check_ref_exist(params.genome, 'genome fasta')
gtf = check_ref_exist(params.gtf, 'gtf file')
bed12 = check_ref_exist(params.bed12, 'bed file')
ref_flat = check_ref_exist(params.ref_flat, 'ref flat file')
star_index = check_ref_exist(params.star_index, 'STAR index')
split_gtf_dir = check_ref_exist(params.split_gtf_dir, 'split gtf directory')
gene_ann = check_ref_exist(params.gene_ann, 'Gene annotation')

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
    file "*.fq.gz" into trimmed_reads
    file "${name}.json" into fastp_json_analysis, fastp_json_report
    file "${name}.html" into fastp_html

    script:
    """
    fastp \\
        --in1 ${reads[0]} \\
        --in2 ${reads[1]} \\
        --out1 ${name}.R1.trimmed.fq.gz \\
        --out2 ${name}.R2.trimmed.fq.gz \\
        --json ${name}.json \\
        --html ${name}.html    
    """
}  


/*
* rm rRNA
*/
process rm_rRNA {

    tag "${name}"

    publishDir "${params.outdir}/${params.proj_name}/data/cleandata/", mode: 'copy',
        saveAs: {filename -> 
            if (filename.indexOf("fq.gz") > 0)  "${filename}"
            else if (filename.indexOf("rrna.log") > 0) "../../logs/rm_rRNA/${filename}"
            else null 
            }    

    cpus = 16

    input:
    file reads from trimmed_reads
    file index from rrna_idx.collect()

    output:
    file "*.fq.gz" into clean_reads, quant_reads
    file "${name}.rrna.log" into rrna_rm_log


    script:
    name = reads[0].toString() - '.R1.trimmed.fq.gz'
    strand_sp = params.bowtie2_lib[params.strand]

    if (params.rm_rrna)
        """
        bowtie2 --reorder --very-sensitive-local ${strand_sp} \\
            -x ${rrna_idx_base} \\
            -1 ${reads[0]} \\
            -2 ${reads[1]} \\
            -p ${task.cpus} \\
            --un-conc-gz ${name}.R%.clean.fq.gz \\
            --no-unal \\
            --no-head \\
            --no-sq \\
            -S ${name}.sam \\
            > ${name}.rrna.log 2>&1
        """
    else
        """
        ln -s ${reads[0]} ${name}.R1.clean.fq.gz
        ln -s ${reads[1]} ${name}.R2.clean.fq.gz

        touch ${name}.rrna.log 2>&1
        """

}


/*
* fq quality control
*/
process reads_qc_summary {

    publishDir "${params.outdir}/${params.proj_name}/result/data_info/", mode: 'copy'

    input:
    file 'fastp_json/*' from fastp_json_analysis.collect()
    file 'rrna/*' from rrna_rm_log.collect()
    
    output:
    file "data_summary.csv" into data_summary, data_summary_qc
    file "reads_gc" into gc_plot, gc_plot_qc
    file "reads_quality" into rq_plot, rq_plot_qc
    file "reads_filter" into rf_plot, rf_plot_qc
    
    script:
    """
    python ${script_dir}/reads_qc/extract_fastp_info.py \\
        --fastp-dir fastp_json \\
        --rrna_dir rrna \\
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
    
    input:
    file reads from clean_reads
    file index from star_index
    
    output:
    file "${sample_name}.bam" into unsort_bam
    file "${sample_name}.Log.final.out" into star_log, star_log_qc
    file "${sample_name}.mapping_rate*" into mapping_plt, mapping_plt_qc

    script:
    sample_name = reads[0].toString() - '.R1.clean.fq.gz'
    if (params.strand == 'unstranded') {
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
    picard FixMateInformation \\
        I=${bam} \\
        O=${name}.sorted.bam \\
        SORT_ORDER=coordinate CREATE_INDEX=true
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
    strand_sp = params.picard_lib[params.strand]
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
    picard CollectInsertSizeMetrics \\
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

    python ${script_dir}/report/report.py .

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
    st_direction = params.stringtie_lib[params.strand]
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

    publishDir "${params.outdir}/${params.proj_name}", mode: 'copy',
            saveAs: {filename ->
                if (filename == "assembly.gtf") "result/assembly/$filename"
                else null
            }

    input:
    file "gtf/*" from assembly_gtf.collect()
    file fasta from genome
    file gtf from gtf
    
    output:
    file "novel.raw.gtf" into novel_gtf, novel_anno_gtf
    file "novel.fa" into novel_fasta
    file "assembly.gtf" into kallisto_idx_gtf, quant_gtf, rmats_gtf
    file "gene_trans.map" into gene2tr_quant, gene2tr_anno
    file "assembly_summary.p*" into assembly_unannotated_fig

    when:
    params.pipeline

    script:
    stranded_flag = params.stranded_flag[params.strand]
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
        --outfile novel.raw.gtf

    gffread novel.raw.gtf \\
        -g ${fasta} \\
        -w novel.fa    

    python ${script_dir}/assembly/get_gene_to_trans.py \\
        --gff assembly.gtf \\
        --out_dir . \\
        --ref
    """
}

/*
* lncRNA prediction
*/ 

process lncRNA_predict {

    publishDir "${params.outdir}/${params.proj_name}/result/lncRNA_prediction/" , mode: 'copy'

    input:
    file novel_fasta from novel_fasta
    file novel_gtf from novel_gtf
    file split_gtf_dir from split_gtf_dir
    file genome from genome
    file gene_ann from gene_ann

    output:
    file "*gene_number.{pdf,png}" into lnc_num_plt
    file "lncRNA_feature*.{pdf,png}" into lnc_feature_plt
    file "lncRNA.gtf" into lncRNA_gtf
    file "TUCP.gtf" into tucp_gtf
    file "Gene_feature.csv" into gene_feature, gf_diff
    file "*fa" into lnc_tucp_fasta

    when:
    params.pipeline

    script:

    """
    /public/python_env/oms_pub/bin/python \\
        /public/software/CPC2/CPC2-beta/bin/CPC2.py \\
            -i ${novel_fasta}

    /public/python_env/work_py3/bin/python ${script_dir}/lncrna/cpc2lnc.py \\
        --cpc cpc2output.txt --gtf ${novel_gtf} \\
        --split_gtf_dir ${split_gtf_dir} \\
        --outfile novel.gtf
    
    /public/python_env/work_py3/bin/python ${script_dir}/lncrna/lnc_feature.py \\
        --novel-lnc novel.gtf \\
        --gtf-split-dir ${split_gtf_dir} \\
        --outdir . \\
        --gene_ann ${gene_ann}

    /public/python_env/work_py3/bin/python ${script_dir}/lncrna/split_gtf_by_type.py \\
        --gtf novel.gtf \\
        --outdir . --novel

    gffread lncRNA.gtf \\
        -g ${genome} \\
        -w lncRNA.fa  

    gffread TUCP.gtf \\
        -g ${genome} \\
        -w TUCP.fa  
    """
}



/*
* quantification
*/
process mk_kallisto_index {
    tag "KALLISTO index"

    publishDir "${params.outdir}/${params.proj_name}", mode: 'copy',
        saveAs: {filename ->
            if (filename == "assembly.fa.kallisto_idx") "data/ref/$filename"
            else null
        }

    input:
    file gtf from kallisto_idx_gtf
    file fasta from genome
    
    output:
    file "assembly.fa" into merged_fa
    file "assembly.fa.kallisto_idx" into kallisto_idx
    
    when:
    params.pipeline

    script:
    """
    gffread ${gtf} \\
        -g ${fasta} \\
        -w assembly.fa

    kallisto index \\
        -i assembly.fa.kallisto_idx \\
        assembly.fa
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
    st_direction = params.kallisto_lib[params.strand]
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
    file gene2tr from gene2tr_anno
    file gene_feature from gene_feature

    output:
    file 'expression_summary' into expression_summary
    file 'deg_input.RData' into deg_obj
    file "expression_summary/*correlation_heatmap.png" into cor_plot
    file "expression_summary/*PCA_plot.png" into pca_plot

    when:
    params.pipeline
    
    script:
    """   
    /public/software/R/R-3.5.1/executable/bin/Rscript ${script_dir}/quant/lnc_kallisto_to_table.R \\
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
            else "${filename}"
        }
        
    input:
    file compare from diff_compare
    file deg_obj from deg_obj
    file sample_group from sample_group
    file gf_diff from gf_diff
    
    output:
    file compare into diff_out_go, diff_out_kegg, diff_out_all, diff_tf, diff_out_ppi, diff_out_enrich_plot, pathway_compare, diff_summary
    file "${compare}/protein_coding.${compare}.Volcano_plot.png" into pcg_diff_plt
    file "${compare}/lncRNA.${compare}.Volcano_plot.png" into lnc_diff_plt
    
    when:
    params.pipeline
    
    script:
    """
    mv ${compare} ${compare}_compare

    /public/software/R/R-3.5.1/executable/bin/Rscript ${script_dir}/quant/diff_edgeR.R \\
        --deg_rdata ${deg_obj} \\
        --compare ${compare} \\
        --sample_inf ${sample_group} \\
        --out_dir ${compare} \\
        --qvalue ${params.exp_diff_pval} \\
        --logfc ${params.exp_lgfc} \\
        --gene_ann ${gf_diff}
    """
}

/*
* diff summary
*/
process diff_exp_summary {

    publishDir "${params.outdir}/${params.proj_name}/result/quantification/expression_summary", mode: 'copy',
        saveAs: {filename -> 
            if (filename.indexOf("png") > 0) null
            else "${filename}"
        }

    input:
    file expression_summary from expression_summary
    file "diff/*" from diff_summary.collect()
    file sample_group from sample_group

    output:
    file "lncRNA" into lnc_diff_summary_out
    file "protein_coding" into pcg_diff_summary_out
    file "lncRNA/*png" into lnc_diff_cls_plt
    file "protein_coding/*png" into pcg_diff_cls_plt

    when:
    params.pipeline

    script:
    """
    /public/software/R/R-3.5.1/executable/bin/Rscript ${script_dir}/quant/quant_report.R \\
        --exp_dir ${expression_summary} \\
        --diff_dir diff \\
        --sample_inf ${sample_group} \\
        --out_dir ./protein_coding \\
        --type protein_coding

    /public/software/R/R-3.5.1/executable/bin/Rscript ${script_dir}/quant/quant_report.R \\
        --exp_dir ${expression_summary} \\
        --diff_dir diff \\
        --sample_inf ${sample_group} \\
        --out_dir ./lncRNA \\
        --type lncRNA
    """
}


/*
* qc report
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
    file ('lnc_number/*') from lnc_num_plt
    file ('lnc_feature/*') from lnc_feature_plt
    file expression_summary from expression_summary
    file ('pca_plot/*') from pca_plot.collect()
    file ('cor_plot/*') from cor_plot.collect()
    file ('pcg_volcano/*') from pcg_diff_plt.collect()
    file ('lnc_volcano/*') from lnc_diff_plt.collect()
    file ('pcg_heatmap/*') from pcg_diff_cls_plt
    file ('lnc_heatmap/*') from lnc_diff_cls_plt

    output:
	file "*.{csv,pdf,png}" into analysis_results
    file "report" into analysis_report

    when:
    params.pipeline

    script:
    stranded_flag = params.stranded_flag[params.strand]
    """
    python ${script_dir}/rnaseq_qc/star_mapping_summary.py \\
        ./star ./

    python ${script_dir}/rnaseq_qc/rnaseq_matrix_summary.py \\
        -d ./rnaseq_matrix -o ./ -m rnaseq

    python ${script_dir}/assembly/assembly_stats.py stats_all \\
        --tmap_dir assembly_tmap ${stranded_flag} \\
        --out_dir .

    python ${script_dir}/report/report.py .

    """
    
}