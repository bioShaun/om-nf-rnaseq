suppressMessages(library(argparser))
suppressMessages(library(plyr))
suppressMessages(library(dplyr))
suppressMessages(library(data.table))
suppressMessages(library(omplotr))

p <- arg_parser("for expression analysis report plot")
p <- add_argument(p, "--exp_dir", help = "expression summary directory")
p <- add_argument(p, "--diff_dir", help = "diff analysis directory")
p <- add_argument(p, "--sample_inf", help = "sample information with sample names and group names")
p <- add_argument(p, "--out_dir", help = "output directory")
p <- add_argument(p, "--qvalue", help = "diff gene qvalue cutoff", default = 0.05)
p <- add_argument(p, "--logfc", help = "diff gene logfc cutoff", default = 1)
argv <- parse_args(p)

DIFF_HEATMAP_GENE = 30000
CUT_TREE_PER = 20
MIN_CLUSTER_NUM = 10
MIN_CLUSTER_POR = 0.005

# read parameters
exp_dir <- argv$exp_dir
diff_dir <- argv$diff_dir
sample_inf <- argv$sample_inf
out_dir <- argv$out_dir
qvalue <- argv$qvalue
logfc <- argv$logfc

samples <- read.delim(sample_inf, stringsAsFactors = F, header = F)
colnames(samples) <- c("sample", "condition")
compare_names <- basename(list.dirs(diff_dir, recursive=F))
diff_df_list = list()
plot_number = min(c(16, length(compare_names)))

diff_genes <- c()
for (each_compare in compare_names) {
  each_compare_name = paste(each_compare, "edgeR.DE_results.txt", sep = ".")
  each_compare_file = file.path(diff_dir, each_compare, each_compare_name)
  each_compare_df <- fread(each_compare_file)
  each_compare_diff_df <- filter(each_compare_df, abs(logFC) >= logfc, FDR <= qvalue)
  each_diff_genes <- each_compare_diff_df$Gene_ID
  diff_genes <- c(diff_genes, each_diff_genes)
}

diff_gene_id <- unique(diff_genes)
raw_gene_count_matrix <- file.path(exp_dir, "Gene.raw.count.txt")
tmm_gene_count_matrix <- file.path(exp_dir, "Gene.TMM.count.txt")
gene_tpm_matrix <- file.path(exp_dir, "Gene.tpm.txt")
raw_gene_count_matrix_df <- fread(raw_gene_count_matrix, check.names = F)
raw_gene_count_matrix_df <- data.frame(raw_gene_count_matrix_df, check.names = F)
tmm_gene_count_matrix_df <- fread(tmm_gene_count_matrix, check.names = F)
tmm_gene_count_matrix_df <- data.frame(tmm_gene_count_matrix_df, check.names = F)
gene_tpm_matrix_df <- fread(gene_tpm_matrix, check.names = F)
gene_tpm_matrix_df <- data.frame(gene_tpm_matrix_df, check.names = F)
out_diff_gene_raw_count_matrix <- filter(raw_gene_count_matrix_df, Gene_ID %in% diff_gene_id)
out_diff_gene_tmm_count_matrix <- filter(tmm_gene_count_matrix_df, Gene_ID %in% diff_gene_id)
out_diff_gene_tpm_matrix <- filter(gene_tpm_matrix_df, Gene_ID %in% diff_gene_id)

if (dim(out_diff_gene_raw_count_matrix)[1] > 0) {
  write.table(out_diff_gene_raw_count_matrix, 
    file = paste(out_dir, "Diff.genes.raw.count.txt",
    sep = "/"), quote = F, row.names = F, sep = "\t")
  write.table(out_diff_gene_tmm_count_matrix, 
    file = paste(out_dir, "Diff.genes.tmm.count.txt",
    sep = "/"), quote = F, row.names = F, sep = "\t")
  write.table(out_diff_gene_tpm_matrix, 
    file = paste(out_dir, "Diff.genes.tpm.txt",
    sep = "/"), quote = F, row.names = F, sep = "\t")
  diff_gene_tpm_matrix <- out_diff_gene_tpm_matrix[, -1]
  rownames(diff_gene_tpm_matrix) <- out_diff_gene_tpm_matrix[, 1]

  # plot heatmap
  diff_gene_count = dim(diff_gene_tpm_matrix)[1]
  if (diff_gene_count <= DIFF_HEATMAP_GENE) {
    om_heatmap(diff_gene_tpm_matrix, samples=samples, outdir=out_dir)
  } else {
    gene_mean_exp <- sort(rowMeans(diff_gene_tpm_matrix), decreasing = T)
    top_genes <- names(gene_mean_exp[1:DIFF_HEATMAP_GENE])
    diff_gene_tpm_matrix <- diff_gene_tpm_matrix[top_genes, ]
    om_heatmap(diff_gene_tpm_matrix, samples, outdir=out_dir)
    diff_gene_count <- dim(diff_gene_tpm_matrix)[1]
  }

  # cluster plot
  cluster_data_dir <- file.path(out_dir, "cluster_data")
  dir.create(cluster_data_dir, showWarnings = F)
  diff_matrix <- as.matrix(diff_gene_tpm_matrix)
  log_diff_matrix <- log2(diff_matrix + 1)
  # center rows, mean substracted
  scale_log_diff_matrix = t(scale(t(log_diff_matrix), scale = F))

  # gene clustering according to centered distance values.
  gene_dist = dist(scale_log_diff_matrix, method = "euclidean")
  hc_genes = hclust(gene_dist, method = "complete")
  gene_partition_assignments <- cutree(as.hclust(hc_genes), h = CUT_TREE_PER/100 * max(hc_genes$height))

  max_cluster_count = max(gene_partition_assignments)
  cluster_num_cutoff = max(c(MIN_CLUSTER_NUM, diff_gene_count * MIN_CLUSTER_POR))

  all_partition_list <- list()
  m = 1
  for (i in 1:max_cluster_count) {
    partition_i = (gene_partition_assignments == i)
    partition_data = scale_log_diff_matrix[partition_i, , drop = F]
    cluster_name <- paste("cluster", i, sep = "_")
    partition_data_df <- as.data.frame(partition_data)
    partition_data_df <- cbind(Gene_id = rownames(partition_data_df), partition_data_df)
    write.table(partition_data_df, file = paste(cluster_data_dir, "/", cluster_name,
      ".txt", sep = ""), quote = F, row.names = F, sep = "\t")
    partition_data_df$cluster <- cluster_name
    melt_partition_data_df <- melt(partition_data_df, id = c("cluster", "Gene_id"))
    out_prefix <- file.path(cluster_data_dir, cluster_name)
    om_cluster_plot(melt_partition_data_df, out_prefix = out_prefix)
    if (dim(partition_data)[1] > cluster_num_cutoff) {
      all_partition_list[[m]] <- partition_data_df
      m <- m + 1
    }
  }

  all_cluster_df <- ldply(all_partition_list, data.frame)
  colnames(all_cluster_df) <- colnames(partition_data_df)
  melt_all_cluster_df <- melt(all_cluster_df, id = c("cluster", "Gene_id"))
  melt_all_cluster_df$variable <- factor(melt_all_cluster_df$variable, levels = samples$sample)
  melt_all_cluster_df$cluster <- factor(melt_all_cluster_df$cluster, levels = unique(melt_all_cluster_df$cluster))
  all_cluster_prefix = file.path(out_dir, 'Diff.genes.cluster')
  om_cluster_plot(melt_all_cluster_df, out_prefix = all_cluster_prefix)
}
