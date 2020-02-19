suppressMessages(library(argparser))
suppressMessages(library(tximport))
suppressMessages(library(tibble))
suppressMessages(library(omplotr))
suppressMessages(library(edgeR))
suppressMessages(library(stringr))
options(stringsAsFactors = F)
options(bitmapType = "cairo")

p <- arg_parser("read kallisto quant files generate expression matrix")
p <- add_argument(p, '--kallisto_dir',  help = 'kallisto quantification directory')
p <- add_argument(p, '--sample_inf', help = 'sample information with sample names and group names')
p <- add_argument(p, '--gene2tr',    help = 'gene id and transcript id mapping file')
p <- add_argument(p, '--out_dir',    help = 'diff analyssi output directory')
p <- add_argument(p, '--gene_feature', help = 'gene feature file.')
argv <- parse_args(p)


om_lnc_boxplot <- function(exp_data, samples, gene_feature,
                           outdir=NULL) {
  
  # normalize data
  plot_data <- norm_exp_data(exp_data)
  
  # reshape matrix and add group inf
  data.m <- reshape2::melt(as.matrix(plot_data))
  colnames(data.m) <- c('gene_id', 'sample_id', 'exp_level')
  row.names(samples) <- samples$sample
  data.m$sample_id <- as.character(data.m$sample_id)
  data.m$group <- samples[data.m$sample_id, 'condition']
  data.m$sample_id <- factor(data.m$sample_id, levels = samples$sample)
  data.m$exp_level <- as.numeric(data.m$exp_level)
  data.m.f <- merge(data.m, gene_feature, by = 'gene_id')
  
  # get sample & group colors
  group_cols <- colorRampPalette(heatmap_col)(length(unique(data.m$group)))
  sample_cols <- colorRampPalette(sample_cols)(length(samples$sample))
  sample_num = length(unique(data.m$sample_id))
  
  theme_set(theme_onmath() +
              theme(axis.text.x = element_text(angle = 90,
                                               vjust = 0.5, hjust = 0),
                    legend.text = element_text(size = rel(0.8)),
                    legend.key = element_blank(),
                    legend.position="top"))
  
  
  boxplot <- ggplot(data.m.f, aes(x = sample_id, y = exp_level, fill = gene_biotype)) +
    geom_boxplot() +
    guides(fill = guide_legend(title = "")) +
    scale_fill_manual(values = group_cols) +
    ylab("Log2(TPM)") + xlab('')



  if (! is.null(outdir)) {
    out_prefix = file.path(outdir,
                           'gene_expression_boxplot')
    plot_width = 4 + sample_num/4
    plot_height = 6 + sample_num/8
    save_ggplot(boxplot, out_prefix,
                width = plot_width,
                height = plot_height)
  } else {
    return(out_plot)
  }
  
}


om_lnc_pca_plot <- function(exp_data, samples, outdir=NULL,
                        out_prefix=NULL, title="") {
  
  plot_data <- norm_exp_data(exp_data)
  PCA_data_mat <- t(apply(plot_data[, 1:dim(plot_data)[2]], 2, as.numeric))
  PCA <- prcomp(PCA_data_mat)
  Summary_PCA <- summary(PCA)
  
  PCA_data <- as.data.frame(PCA$x[, 1:dim(samples)[1]])
  sample_name <- rownames(PCA$x)
  match_index <- match(sample_name, samples$sample, nomatch = 0)
  group_name <- samples$condition[match_index]
  PCA_data$Sample <- sample_name
  PCA_data$Group <- group_name
  
  sample_types <- as.character(unique(samples$condition))
  nsamples <- length(sample_types)
  sample_colors <- colorRampPalette(heatmap_col)(nsamples)
  
  pca_plot <- ggplot(PCA_data, aes(PC1, PC2)) +
    geom_point(aes(colour = Group), size = rel(3)) +
    geom_text(aes(label = Sample), vjust = 0,
              hjust = 0.5, color = "black",
              size = rel(2)) +
    geom_vline(xintercept = 0, linetype = 2,
               colour = "grey60", size = rel(0.5)) +
    geom_hline(yintercept = 0, linetype = 2,
               colour = "grey60", size = rel(0.5)) +
    theme_onmath() +
    scale_color_manual(values = sample_colors) +
    labs(title = title,
         x = paste0("PC1: ", Summary_PCA$importance[2, 1] * 100, "% variance"),
         y = paste0("PC2: ", Summary_PCA$importance[2, 2] * 100, "% variance"))
  
  if (! is.null(outdir)) {
    out_prefix <- file.path(outdir, "PCA_plot")
    save_ggplot(pca_plot, out_prefix)
  } else if (! is.null(out_prefix)) {
    out_dir <- dirname(out_prefix)
    path_name <- save_mkdir(out_dir)
    save_ggplot(ggplot_out = pca_plot,
                output = out_prefix)
  }
  
  return(pca_plot)
}


om_lnc_correlation_plot <- function(exp_data, samples, 
                                    outdir=NULL, title="") {
  data <- norm_exp_data(exp_data)
  sample_types <- as.character(unique(samples$condition))
  rep_names <- as.character(samples$sample)
  data <- data[, colnames(data) %in% samples$sample, drop = F]
  nsamples <- length(sample_types)
  sample_colors <- colorRampPalette(heatmap_col)(nsamples)
  data <- as.matrix(data)
  sample_cor <- cor(data, method = "pearson", use = "pairwise.complete.obs")
  sample_cor_df <- as.data.frame(sample_cor)
  sample_cor_df <- cbind(Sample = rownames(sample_cor_df), sample_cor_df)
  
  Group <- sample_colors[1:length(unique(samples$condition))]
  names(Group) <- unique(samples$condition)
  ann_color = data.frame(group = samples$condition)
  rownames(ann_color) <- samples$sample
  ann_colors = list(group = Group)
  
  theme_set(theme_onmath() + theme(legend.position = c(0.5, 0.5)))
  sample_num = length(colnames(exp_data))
  heatmap_width <- (sample_num - 5)/5 + 8
  heatmap_heigh <- (sample_num - 5)/5 + 6
  fontsize = (sample_num - 5)/10 + 7
  cellwidth <- (heatmap_width - 1) * 30/sample_num
  
  draw_heatmap <- function(){
    pheatmap::  pheatmap(sample_cor, annotation_col = ann_color,
                         annotation_colors = ann_colors,
                         annotation_row = ann_color, annotation_names_row = F,
                         annotation_names_col = F,
                         color = rev(cor_plot_col), border_color = NA,
                         cellwidth = cellwidth, fontsize = fontsize,
                         cellheight = cellwidth,
                         main=title)
  }
  
  if (! is.null(outdir)) {
    
    out_prefix <- file.path(outdir,
                            paste(title, 'sample_correlation_heatmap', sep='_'))
    write.table(sample_cor_df, file = paste(out_prefix, "txt",sep = "."),
                quote = F, sep = "\t", row.names = F)
    save_general_plot(draw_heatmap(), out_prefix,
                      width = heatmap_width,
                      height = heatmap_heigh,
                      plot_type='pdf')
    save_general_plot(draw_heatmap(), out_prefix,
                      width = heatmap_width,
                      height = heatmap_heigh,
                      plot_type='png')
    
  }
  
  return(draw_heatmap())
  
}

## read parameters
sample_inf <- argv$sample_inf
kallisto_dir <- argv$kallisto_dir
gene2tr_file <- argv$gene2tr
expression_stat_dir <- argv$out_dir
gene_feature <- argv$gene_feature

## directory prepare
dir.create(expression_stat_dir, showWarnings = FALSE)

samples <- read.delim(sample_inf, stringsAsFactors = F, header = F)
colnames(samples) <- c('sample', 'condition')
files <- file.path(kallisto_dir, samples$sample, "abundance.h5")
names(files) <- samples$sample
gene2tr <- read.delim(gene2tr_file, header = FALSE)
colnames(gene2tr) <- c('gene_id', 'transcript_id')
tx2gene <- gene2tr[,c('transcript_id', 'gene_id')]


## normalized expression matrix
txi <- tximport(files, type = "kallisto", tx2gene = tx2gene)
cts <- txi$counts
y <- DGEList(cts)
normfactors <- calcNormFactors(y)
norm_cts <- sweep(cts, 2, normfactors$samples$norm.factors, FUN = '/')
gene_tpm_matrix <- txi$abundance

normMat <- txi$length
normMat <- normMat/exp(rowMeans(log(normMat)))
o <- log(calcNormFactors(cts/normMat)) + log(colSums(cts/normMat))
deg_obj <- DGEList(cts, group=samples$condition)
deg_obj <- scaleOffset(deg_obj, t(t(log(normMat)) + o))
#keep <- filterByExpr(deg_obj)
# bug: In min(n[n > 1L]) : no non-missing arguments to min; returning Inf
keep <- rowSums(cpm(y) > 0.1) >= 1 
deg_obj_list <- list(deg_obj=deg_obj[keep, ],
                     normfactors=normfactors$samples$norm.factors)
save(deg_obj_list, file='deg_input.RData')

## output quant table
output_matrix <- function(exp_matrix, gene_ann, out_dir, out_name) {
    out_df <- as.data.frame(exp_matrix)
    out_df <- round(out_df, 3)
    out_df <- rownames_to_column(out_df, var="Gene_ID")
    out_df <- merge(gene_ann, out_df, by.x='gene_id', by.y='Gene_ID', all.y=T)
    colnames(out_df)[1] <- 'Gene_ID'
    out_df[is.na(out_df)] <- '--'
    write.table(out_df,
                file = paste(out_dir, out_name, sep = '/'),
                quote=F, row.names = F, sep = '\t')
}

## boxplot
gf_df <- read.csv(gene_feature)
gf_df <- subset(gf_df, select = -c(Isoform_number, Length, Exon_number))

output_matrix(cts, gf_df, expression_stat_dir, 'Gene.raw.count.txt')
output_matrix(gene_tpm_matrix, gf_df, expression_stat_dir, 'Gene.tpm.txt')
output_matrix(norm_cts, gf_df, expression_stat_dir, 'Gene.TMM.count.txt')

gf_df <- gf_df[str_detect(gf_df$gene_biotype, 'protein_coding') | str_detect(gf_df$gene_biotype, 'lncRNA'), ]
gf_df[str_detect(gf_df$gene_biotype, 'lncRNA'),]$gene_biotype <- 'lncRNA'
om_lnc_boxplot(gene_tpm_matrix, samples, gf_df,
           outdir=expression_stat_dir)

## PCA plot & sample correlation

for (gt in unique(gf_df$gene_biotype)) {
  gt_genes_df <- dplyr::filter(gf_df, gene_biotype == gt)
  gt_genes <- as.vector(gt_genes_df$gene_id)
  gt_gene_tpm_matrix <- gene_tpm_matrix[rownames(gene_tpm_matrix) %in% gt_genes,]
  out_prefix <- file.path(expression_stat_dir,
                          paste(gt, 'PCA_plot', sep = '_'))
  om_lnc_pca_plot(gt_gene_tpm_matrix, samples, out_prefix = out_prefix, title = gt)
  om_lnc_correlation_plot(gt_gene_tpm_matrix, samples, outdir = expression_stat_dir, title = gt)
}


