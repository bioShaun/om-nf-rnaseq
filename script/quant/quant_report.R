suppressMessages(library(stringr))
suppressMessages(library(argparser))
suppressMessages(library(plyr))
suppressMessages(library(dplyr))
suppressMessages(library(data.table))
suppressMessages(library(omplotr))

options(stringsAsFactors = F)

p <- arg_parser("for expression analysis report plot")
p <- add_argument(p, "--exp_dir", help = "expression summary directory")
p <- add_argument(p, "--diff_dir", help = "diff analysis directory")
p <- add_argument(p, "--sample_inf", help = "sample information with sample names and group names")
p <- add_argument(p, "--out_dir", help = "output directory")
p <- add_argument(p, "--qvalue", help = "diff gene qvalue cutoff", default = 0.05)
p <- add_argument(p, "--logfc", help = "diff gene logfc cutoff", default = 1)
p <- add_argument(p, "--grp_num", help = "minimal group number to perform cluster analysis", default = 3)
argv <- parse_args(p)


	om_lnc_heatmap <- function(exp_data, samples,
			       outdir=NULL, out_prefix=NULL,
			       scale='row',
			       cluster_rows=T,
			       cluster_cols=T,
			       title="") {

	  # normalize exp data
	  plot_data <- norm_exp_data(exp_data)

	  # heatmap outer bar color
	  if (length(unique(samples$condition)) > length(heatmap_col)) {
	    heatmap_col <- colorRampPalette(heatmap_col)(length(unique(samples$condition)))
	  }

	  # heatmap cell color
	  cell_cols <- RColorBrewer::brewer.pal(10,"RdYlGn")
	  gradiant_cell_cols <- rev(colorRampPalette(cell_cols)(100))

	  # map group to bar color
	  Group <- heatmap_col[1:length(unique(samples$condition))]
	  names(Group) <- unique(samples$condition)
	  ann_color = data.frame(group = samples$condition)
	  rownames(ann_color) <- samples$sample
	  ann_colors = list(group = Group)

	  # adjust output height & width
	  sample_num = length(colnames(plot_data))
	  heatmap_width <- (sample_num - 5)/3 + 2.5
	  heatmap_heigh <- (sample_num - 5)/3 + 4.5
	  fontsize = (sample_num - 5)/10 + 4.5
	  cellwidth <- (heatmap_width - 0.5) * 50/sample_num

	  ## TODO legend size
	  draw_heatmap <- function(){
	    pheatmap::pheatmap(plot_data, show_rownames = F,
			       annotation_col = ann_color, annotation_colors = ann_colors,
			       annotation_legend = T, annotation_names_col = F,
			       color = gradiant_cell_cols,
			       treeheight_row = 0, scale = scale, fontsize = fontsize,
			       cellwidth = cellwidth, border_color = NA,
			       cluster_rows=cluster_rows,
			       cluster_cols=cluster_cols,
			       main=title)
	  }

	  if (! is.null(outdir)) {
	    out_prefix <- file.path(outdir, 'Diff_genes_heatmap')
	    save_general_plot(draw_heatmap(), out_prefix,
			      width = heatmap_width,
			      height = heatmap_heigh,
			      plot_type='pdf')
	    save_general_plot(draw_heatmap(), out_prefix,
			      width = heatmap_width,
			      height = heatmap_heigh,
			      plot_type='png')
	  } else if (! is.null(out_prefix)) {
	    out_dir <- dirname(out_prefix)
	    path_name <- save_mkdir(out_dir)
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


	gene_type_heatmap <- function (gene_type, exp_df, samples, out_dir) {
	  out_diff_gene_tpm_matrix <- exp_df[str_detect(exp_df$gene_biotype, gene_type),]
	  diff_gene_tpm_matrix <- out_diff_gene_tpm_matrix[, -1]
	  rownames(diff_gene_tpm_matrix) <- out_diff_gene_tpm_matrix[, 1]
	  diff_gene_tpm_matrix <- diff_gene_tpm_matrix[, samples$sample]

	  # plot heatmap
	  heatmap_name = paste(gene_type, "diff_genes_heatmap", sep='_')
	  heatmap_prefix = file.path(out_dir, heatmap_name)
	  diff_gene_count = dim(diff_gene_tpm_matrix)[1]
	  if (diff_gene_count <= DIFF_HEATMAP_GENE) {
	    om_lnc_heatmap(diff_gene_tpm_matrix, samples=samples, out_prefix=heatmap_prefix, title=gene_type)
	  } else {
	    gene_mean_exp <- sort(rowMeans(diff_gene_tpm_matrix), decreasing = T)
	    top_genes <- names(gene_mean_exp[1:DIFF_HEATMAP_GENE])
	    diff_gene_tpm_matrix <- diff_gene_tpm_matrix[top_genes, ]
	    om_lnc_heatmap(diff_gene_tpm_matrix, samples, out_prefix=heatmap_prefix, title=gene_type)
	    diff_gene_count <- dim(diff_gene_tpm_matrix)[1]
	  }

	}


	om_lnc_cluster_plot <- function(plot_data, out_prefix=NULL) {

	  #cluster_number = length(unique(plot_data$gene_biotype))
	  #col_theme <- colorRampPalette(heatmap_col)(cluster_number)
	  theme_cluster <- theme_onmath() + theme(axis.text.x = element_text(vjust = -0.2,
									     angle = 90,
									     size = rel(0.8)))

	  cluster_plot <- ggplot(plot_data, aes(x=variable, y=value, group = Gene_ID, color=gene_biotype)) +
	    geom_line(alpha = 0.5) +
	    scale_color_brewer(palette = 'Set1') +
	    xlab("") + ylab("Scaled log10(tpm+1)") + theme_cluster

	  if (! is.null(out_prefix)) {
	    plot_height <- 6  
	    plot_width <- 8 
	    save_ggplot(cluster_plot, out_prefix,
			width = plot_width, height = plot_height)
	  }

	  return(cluster_plot)

	}


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
	min_grp_num <- argv$grp_num

	dir.create(out_dir, showWarnings = FALSE)

	samples <- read.delim(sample_inf, stringsAsFactors = F, header = F)
	colnames(samples) <- c("sample", "condition")
	compare_names <- basename(list.dirs(diff_dir, recursive=F))
	diff_df_list = list()
	plot_number = min(c(16, length(compare_names)))
	grp_number <- length(unique(samples$condition))


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
	    file = paste(out_dir, "/", "Diff_genes_raw_count.txt",
	    sep = ""), quote = F, row.names = F, sep = "\t")
	  write.table(out_diff_gene_tmm_count_matrix, 
	    file = paste(out_dir, "/", "Diff_genes_tmm_count.txt",
	    sep = ""), quote = F, row.names = F, sep = "\t")
	  write.table(out_diff_gene_tpm_matrix, 
	    file = paste(out_dir, "/", "Diff_genes_tpm.txt",
	    sep = ""), quote = F, row.names = F, sep = "\t")

	  gene_type_heatmap('protein_coding', out_diff_gene_tpm_matrix, samples, out_dir)
	  gene_type_heatmap('lncRNA', out_diff_gene_tpm_matrix, samples, out_dir)

	  # cluster plot
	  if (grp_number >= min_grp_num) {
	    exp_df <- out_diff_gene_tpm_matrix
	    exp_df <- exp_df[str_detect(exp_df$gene_biotype, 'protein_coding') | str_detect(exp_df$gene_biotype, 'lncRNA'), ]
	    exp_df[str_detect(exp_df$gene_biotype, 'lncRNA'),]$gene_biotype <- 'lncRNA'
	    rownames(exp_df) <- exp_df[, 1]
	    diff_gene_tpm_matrix <- exp_df[, samples$sample]
	    grp_exp_df <- merge(t(diff_gene_tpm_matrix), samples, by.x=0, by.y='sample')
	    grp_exp_df <- grp_exp_df[, -1]
	    grp_mean_df <- aggregate( .~ condition, grp_exp_df, mean)
	    rownames(grp_mean_df) <- grp_mean_df$condition
	    grp_mean_df <- t(grp_mean_df[, -1])
	  
	    cluster_data_dir <- file.path(out_dir, "cluster_plot")
	    dir.create(cluster_data_dir, showWarnings = F)
	    diff_matrix <- as.matrix(grp_mean_df)
	    log_diff_matrix <- log2(diff_matrix + 1)
	    # center rows, mean substracted
	    scale_log_diff_matrix = t(scale(t(log_diff_matrix), scale = F))
	  
	    # gene clustering according to centered distance values.
	    gene_dist = dist(scale_log_diff_matrix, method = "euclidean")
	    hc_genes = hclust(gene_dist, method = "complete")
	    gene_partition_assignments <- cutree(as.hclust(hc_genes), h = CUT_TREE_PER/100 * max(hc_genes$height))
	  
	    max_cluster_count = max(gene_partition_assignments)
	    #cluster_num_cutoff = max(c(MIN_CLUSTER_NUM, diff_gene_count * MIN_CLUSTER_POR))
	  
	    all_partition_list <- list()
	    m = 1
	    for (i in 1:max_cluster_count) {
	      partition_i = (gene_partition_assignments == i)
	      partition_data = scale_log_diff_matrix[partition_i, , drop = F]
	      partition_genes = names(gene_partition_assignments[partition_i])
	      partition_genes_df = filter(exp_df, Gene_ID %in% partition_genes)
	      ann_cols <- colnames(partition_genes_df)[! colnames(partition_genes_df) %in% samples$sample]
	      partition_genes_df = partition_genes_df[, ann_cols]
	  
	      cluster_name <- paste("cluster", i, sep = "_")
	      partition_data_df <- as.data.frame(partition_data)
	      partition_data_df <- cbind(Gene_ID = rownames(partition_data_df), partition_data_df)
	      #write.table(partition_data_df, file = paste(cluster_data_dir, "/", cluster_name,
	      #  ".txt", sep = ""), quote = F, row.names = F, sep = "\t")
	      #partition_data_df$cluster <- cluster_name
	      partition_genes_df$cluster <- cluster_name
	      partition_genes_df <- merge(partition_genes_df, partition_data_df, by="Gene_ID")
	      melt_partition_data_df <- reshape2::melt(partition_genes_df, 
					     id.vars = c("Gene_ID", "gene_biotype"),
					     measure.vars=unique(samples$condition))
	      out_prefix <- file.path(cluster_data_dir, cluster_name)
	      om_lnc_cluster_plot(melt_partition_data_df, out_prefix = out_prefix)
	      all_partition_list[[m]] <- partition_genes_df
	      m <- m + 1
	    }
	  
	    all_cluster_df <- ldply(all_partition_list, data.frame)
	    colnames(all_cluster_df) <- colnames(partition_genes_df)
	    #melt_all_cluster_df <- melt(all_cluster_df, id = c("cluster", "Gene_id"))
	    #melt_all_cluster_df$variable <- factor(melt_all_cluster_df$variable, levels = samples$sample)
	    #melt_all_cluster_df$cluster <- factor(melt_all_cluster_df$cluster, levels = unique(melt_all_cluster_df$cluster))
	    all_cluster_prefix = file.path(out_dir, 'Diff_genes_cluster.txt')
	    write.table(all_cluster_df, file=all_cluster_prefix, quote = F, row.names = F, sep = "\t")
	    #om_cluster_plot(melt_all_cluster_df, out_prefix = all_cluster_prefix)
	  } 
	}
