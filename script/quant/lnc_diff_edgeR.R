suppressMessages(library(argparser))
suppressMessages(library(stringr))
suppressMessages(library(tximport))
suppressMessages(library(edgeR))
suppressMessages(library(dplyr))
suppressMessages(library(rhdf5))
suppressMessages(library(omplotr))
options(stringsAsFactors = F)


p <- arg_parser("perform differential analysis")
p <- add_argument(p, '--deg_rdata', help = 'edgeR DEGList RData')
p <- add_argument(p, '--compare', help = 'compare name')
p <- add_argument(p, '--sample_inf', help = 'sample information with sample names and group names')
p <- add_argument(p, '--out_dir', help = 'output directory')
p <- add_argument(p, '--qvalue', help = 'diff gene qvalue cutoff', default = 0.05)
p <- add_argument(p, '--logfc', help = 'diff gene logfc cutoff',  default = 1)
p <- add_argument(p, '--dispersion', help = 'for analysis with no replicates',  default = 0.1)
p <- add_argument(p, '--gene_ann', help = 'gene annotation')
argv <- parse_args(p)


om_lnc_volcano_plot <- function(dpa_results, compare_name,
                            logfc=1, qvalue=0.05, outdir=NULL,
                            gt="ALL") {
  
  plot_cols <- c('#E41A1C', '#FF7F00', '#984EA3',
                 '#377EB8', '#4DAF4A', '#E5C494',
                 '#999999')
  plot_order <- c('Up-regulated-protein_coding',
                        'Up-regulated-lncRNA',
                        'Up-regulated-others',
                        'Down-regulated-protein_coding',
                        'Down-regulated-lncRNA',
                        'Down-regulated-others',
                        'Non-significant')
  names(plot_cols) <- plot_order
  dpa_results$logFDR <- -log10(dpa_results$FDR)
  dpa_results$status <- "Non-significant"
  dpa_results[dpa_results$FDR <= qvalue & dpa_results$logFC >= logfc, ]$status <- 'Up-regulated'
  dpa_results[dpa_results$FDR <= qvalue & dpa_results$logFC <= -logfc, ]$status <- 'Down-regulated'
  dpa_results$merge_status <- paste(
    dpa_results$status, 
    dpa_results$gene_biotype, 
    sep='-')
  dpa_results[str_detect(dpa_results$merge_status, 'Non-significant'), ]$merge_status <- 'Non-significant'
  diff_count <- table(dpa_results$merge_status)
  diff_label <- paste(names(diff_count), diff_count, sep=': ')
  names(diff_label) <- names(diff_count)
  dpa_results$diff_num <- diff_count[dpa_results$merge_status]
  dpa_results$plot_label <- paste(dpa_results$merge_status, dpa_results$diff_num, sep=': ')

  # labels and titles

  logFC_max = max(dpa_results$logFC)
  logFC_min = min(dpa_results$logFC)

  logFC_limit <- ceiling(max(c(abs(logFC_min), logFC_max)))
  logFC_limit <- ifelse(logFC_limit < 8, 8, logFC_limit)
  logFC_limit <- ifelse(logFC_limit > 15, 15, logFC_limit)
  logFDR_limit <- ceiling(max(dpa_results$logFDR))
  logFDR_limit <- ifelse(logFDR_limit > 50, 50, logFDR_limit)
  logFDR_limit <- ifelse(logFDR_limit < 35, 35, logFDR_limit)
  
  theme_set(theme_onmath() +
              theme(legend.key = element_blank(),
                    panel.grid.major = element_blank(),
                    legend.direction = "vertical"))
  
  

  y_line_pos = round(-log10(qvalue), 1)
  dpa_results$merge_status <- factor(dpa_results$merge_status,
                                     levels = plot_order)
  p <- ggplot(dpa_results, aes(logFC, logFDR, colour = merge_status)) +
    geom_point(size=0.6) +
    geom_hline(yintercept = y_line_pos, lty = 4, size = 0.45) +
    geom_vline(xintercept = -(logfc),lty = 4, size = 0.45) +
    geom_vline(xintercept = logfc, lty = 4, size = 0.45) +
    xlab("logFC") + ylab("-log10(FDR)")
  
  p <- p + guides(colour = guide_legend(
    override.aes = list(size=2),
    title = "")) +
    ggtitle(compare_name) + scale_y_continuous(breaks = c(0, 10, 20, 30),
                                               limits = c(0, logFDR_limit)) +
    scale_x_continuous(breaks = c(-8,-4, -2, -1, 0, 1, 2, 4, 8),
                       limits = c(-logFC_limit, logFC_limit)) +
    scale_color_manual(values = plot_cols, labels = diff_label)
    
    
  plot_height <- 6 
  plot_width <- 8 
  if (! is.null(outdir)) {
    if (gt == "ALL") {
       plt_name <- paste(compare_name, 'Volcano_plot', sep = '_')
    } else {
       plt_name <- paste(gt, compare_name, 'Volcano_plot', sep = '_')
    }
    out_prefix <- file.path(outdir, plt_name)
    save_ggplot(p, out_prefix, width = plot_width, height = plot_height)
  }
  return(p)
  
}



## read aguments
deg_rdata <- argv$deg_rdata
compare <- argv$compare
sample_inf <- argv$sample_inf
outdir <- argv$out_dir
qvalue <- argv$qvalue
logfc <- argv$logfc
dispersion <- argv$dispersion
gene_ann <- argv$gene_ann

gf_df <- read.csv(gene_ann)
gf_df <- subset(gf_df, select = -c(Isoform_number, Length, Exon_number))

dir.create(outdir, showWarnings = FALSE)
samples <- read.delim(sample_inf, stringsAsFactors = F, header = F)
colnames(samples) <- c('sample', 'condition')

load(deg_rdata)

each_pair <- unlist(strsplit(compare, split = "_vs_"))
con1_sample <- samples[samples$condition == each_pair[1],'sample']
con2_sample <- samples[samples$condition == each_pair[2],'sample']
each_pair_samples <- samples[samples$sample %in% c(con1_sample, con2_sample),]

## diff
deg_obj <- deg_obj_list$deg_obj
if (length(con1_sample) > 1 && length(con2_sample) > 1) {
  y <- estimateDisp(deg_obj)
  et <- exactTest(y, pair=rev(each_pair))
} else {
  et <- exactTest(deg_obj, pair=rev(each_pair), dispersion=dispersion)
}

tTags <- topTags(et,n=NULL)
new_tTags <- tTags$table
new_tTags <- new_tTags[, !(names(new_tTags) %in% c("logCPM"))]

norm_counts <- sweep(deg_obj$counts, 2, deg_obj_list$normfactors, FUN = '/')
pair_norm_counts = norm_counts[, c(con1_sample, con2_sample)]

merged_df <- merge(pair_norm_counts, new_tTags, by.x = 0, by.y = 0, all.y = T)
colnames(merged_df)[1] <- 'gene_id'
merged_df <- merge(gf_df, merged_df, by='gene_id', all.y=T)
merged_df[is.na(merged_df)] <- '--'
sorted_merged_df <- arrange(merged_df, FDR)
colnames(sorted_merged_df)[1] <- 'Gene_ID'

out_file_name_prefix <- paste(outdir, '/', compare, sep = '')
up_regulate_name_prefix <- paste(outdir, '/', compare, '.','UP',  sep = '')
down_regulate_name_prefix <- paste(outdir, '/', compare, '.','DOWN', sep = '')
pcg_out_file_name_prefix <- paste(outdir, '/', 'protein_coding.', compare, sep = '')
pcg_up_regulate_name_prefix <- paste(outdir, '/', 'protein_coding.', compare, '.','UP',  sep = '')
pcg_down_regulate_name_prefix <- paste(outdir, '/', 'protein_coding.', compare, '.','DOWN', sep = '')
lnc_out_file_name_prefix <- paste(outdir, '/', 'lncRNA.', compare, sep = '')
lnc_up_regulate_name_prefix <- paste(outdir, '/', 'lncRNA.', compare, '.','UP',  sep = '')
lnc_down_regulate_name_prefix <- paste(outdir, '/', 'lncRNA.', compare, '.','DOWN', sep = '')

diff_genes <- c()
up_regulate_df <- filter(sorted_merged_df, logFC >= logfc, FDR <= qvalue)
down_regulate_df <- filter(sorted_merged_df, logFC <= -(logfc), FDR <= qvalue)
diff_genes <- c(diff_genes, up_regulate_df$Gene_ID, down_regulate_df$Gene_ID)
write.table(sorted_merged_df, file=paste(out_file_name_prefix, 'edgeR.DE_results.txt', sep = '.'), sep='\t', quote=F, row.names=F)
lnc_genes <- sorted_merged_df[str_detect(sorted_merged_df$gene_biotype, 'lncRNA'),]$Gene_ID
pcg_genes <- sorted_merged_df[str_detect(sorted_merged_df$gene_biotype, 'protein_coding'),]$Gene_ID

if (dim(up_regulate_df)[1] > 0) {
  write.table(up_regulate_df, file=paste(up_regulate_name_prefix, 'edgeR.DE_results.txt', sep = '.'), sep='\t', quote=F, row.names=F)
  write(as.character(up_regulate_df$Gene_ID), file = paste(up_regulate_name_prefix, 'edgeR.DE_results.diffgenes.txt', sep = '.'), sep = '\n')
  up_pcg_genes <- up_regulate_df$Gene_ID[up_regulate_df$Gene_ID %in% pcg_genes]
  up_lnc_genes <- up_regulate_df$Gene_ID[up_regulate_df$Gene_ID %in% lnc_genes]
  write(as.character(up_pcg_genes), file = paste(pcg_up_regulate_name_prefix, 'edgeR.DE_results.diffgenes.txt', sep = '.'), sep = '\n')
  write(as.character(up_lnc_genes), file = paste(lnc_up_regulate_name_prefix, 'edgeR.DE_results.diffgenes.txt', sep = '.'), sep = '\n')
}
if (dim(down_regulate_df)[1] > 0) {
  write.table(down_regulate_df, file=paste(down_regulate_name_prefix, 'edgeR.DE_results.txt', sep = '.'), sep='\t', quote=F, row.names=F)
  write(as.character(down_regulate_df$Gene_ID), file = paste(down_regulate_name_prefix, 'edgeR.DE_results.diffgenes.txt', sep = '.'), sep = '\n')
  down_pcg_genes <- up_regulate_df$Gene_ID[down_regulate_df$Gene_ID %in% pcg_genes]
  down_lnc_genes <- up_regulate_df$Gene_ID[down_regulate_df$Gene_ID %in% lnc_genes]
  write(as.character(down_pcg_genes), file = paste(pcg_down_regulate_name_prefix, 'edgeR.DE_results.diffgenes.txt', sep = '.'), sep = '\n')
  write(as.character(down_lnc_genes), file = paste(lnc_down_regulate_name_prefix, 'edgeR.DE_results.diffgenes.txt', sep = '.'), sep = '\n')
}

## write diff gene list
if (length(diff_genes) > 0) {
  write(as.character(diff_genes), file = paste(out_file_name_prefix, 'ALL.edgeR.DE_results.diffgenes.txt', sep = '.'), sep = '\n')
  pcg_diff_genes <- diff_genes[diff_genes %in% pcg_genes]
  lnc_diff_genes <- diff_genes[diff_genes %in% lnc_genes]
  write(as.character(pcg_diff_genes), file = paste(pcg_out_file_name_prefix, 'ALL.edgeR.DE_results.diffgenes.txt', sep = '.'), sep = '\n')
  write(as.character(lnc_diff_genes), file = paste(lnc_out_file_name_prefix, 'ALL.edgeR.DE_results.diffgenes.txt', sep = '.'), sep = '\n')
}
## volcano plot
sorted_merged_df[str_detect(sorted_merged_df$gene_biotype, 'lncRNA'),]$gene_biotype <- 'lncRNA'
sorted_merged_df[sorted_merged_df$gene_biotype != 'lncRNA' & sorted_merged_df$gene_biotype != 'protein_coding',]$gene_biotype <- 'others'

om_lnc_volcano_plot(sorted_merged_df, compare, logfc, qvalue, outdir)
pcg_dpa_results <- dplyr::filter(sorted_merged_df, gene_biotype == 'protein_coding')
om_lnc_volcano_plot(pcg_dpa_results, compare, logfc, qvalue, outdir, gt='protein_coding')
lnc_dpa_results <- dplyr::filter(sorted_merged_df, gene_biotype == 'lncRNA')
om_lnc_volcano_plot(lnc_dpa_results, compare, logfc, qvalue, outdir, gt='lncRNA')

