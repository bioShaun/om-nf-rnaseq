suppressMessages(library(argparser))
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
argv <- parse_args(p)


## read aguments
deg_rdata <- argv$deg_rdata
compare <- argv$compare
sample_inf <- argv$sample_inf
outdir <- argv$out_dir
qvalue <- argv$qvalue
logfc <- argv$logfc
dispersion <- argv$dispersion

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
sorted_merged_df <- arrange(merged_df, FDR)
colnames(sorted_merged_df)[1] <- 'Gene_ID'

out_file_name_prefix <- paste(outdir, '/', compare, sep = '')
up_regulate_name_prefix <- paste(outdir, '/', compare, '.',each_pair[1], '-UP',  sep = '')
down_regulate_name_prefix <- paste(outdir, '/', compare, '.',each_pair[2], '-UP', sep = '')

diff_genes <- c()
up_regulate_df <- filter(sorted_merged_df, logFC >= logfc, FDR <= qvalue)
down_regulate_df <- filter(sorted_merged_df, logFC <= -(logfc), FDR <= qvalue)
diff_genes <- c(diff_genes, up_regulate_df$Gene_ID, down_regulate_df$Gene_ID)
write.table(sorted_merged_df, file=paste(out_file_name_prefix, 'edgeR.DE_results.txt', sep = '.'), sep='\t', quote=F, row.names=F)
if (dim(up_regulate_df)[1] > 0) {
  write.table(up_regulate_df, file=paste(up_regulate_name_prefix, 'edgeR.DE_results.txt', sep = '.'), sep='\t', quote=F, row.names=F)
  write(as.character(up_regulate_df$Gene_ID), file = paste(up_regulate_name_prefix, 'edgeR.DE_results.diffgenes.txt', sep = '.'), sep = '\n')
}
if (dim(down_regulate_df)[1] > 0) {
  write.table(down_regulate_df, file=paste(down_regulate_name_prefix, 'edgeR.DE_results.txt', sep = '.'), sep='\t', quote=F, row.names=F)
  write(as.character(down_regulate_df$Gene_ID), file = paste(down_regulate_name_prefix, 'edgeR.DE_results.diffgenes.txt', sep = '.'), sep = '\n')
}
## write diff gene list
if (length(diff_genes) > 0) {
  write(as.character(diff_genes), file = paste(out_file_name_prefix, 'ALL.edgeR.DE_results.diffgenes.txt', sep = '.'), sep = '\n')
}
## volcano plot
om_volcano_plot(sorted_merged_df, compare, logfc, qvalue, outdir)
