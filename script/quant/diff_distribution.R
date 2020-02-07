suppressMessages(library(ggplot2))
suppressMessages(library(dplyr))
suppressMessages(library(plyr))
suppressMessages(library(RColorBrewer))
suppressMessages(library(scales))
suppressMessages(library(argparser))
suppressMessages(library(ggnewscale))
suppressMessages(library(stringr))

options(stringsAsFactors = F)

lnc_distribution_heatmap <- function(diff_df, chrom_df, out_prefix, col_pal) {
  if (dim(diff_df)[1] == 0) {
    return(1)
  }
  plt_title <- basename(out_prefix)
  plt_title <- str_remove(plt_title, '_diff_gene_location')
  max_count <- max(diff_df$gene_count)
  up_df <- filter(diff_df, thickEnd == 'UP')
  down_df <- filter(diff_df, thickEnd == 'DOWN')
  up_df$ymin <- 0
  up_df$ymax <- 1
  down_df$ymin <- -1
  down_df$ymax <- 0
  
  if (col_pal == 1) {
    up_col <- colorRampPalette(brewer.pal(9, 'Reds'))(100)
    down_col <- colorRampPalette(brewer.pal(9, 'Blues'))(100)    
  } else {
    up_col <- colorRampPalette(brewer.pal(9, 'Oranges'))(100)
    down_col <- colorRampPalette(brewer.pal(9, 'Greens'))(100)       
  }
  
  deci_num <- floor(log10(max(chrom_df$end)))
  max_chr_len = ceiling(max(chrom_df$end)/(10^deci_num))
  leading_zero <- paste(rep("0", deci_num - 6), collapse = "")
  
  chrom_order <- chrom_df$chrom
  down_df$chrom <- factor(down_df$chrom,
                          levels = chrom_order)
  up_df$chrom <- factor(up_df$chrom,
                        levels = chrom_order)
  chrom_df$chrom <- factor(chrom_df$chrom,
                           levels = chrom_order)
  
  p <- ggplot() + 
    geom_rect(data = down_df,
              aes(xmin=start, xmax=end, ymin=ymin, 
                  ymax=ymax, fill = gene_count)) +
    scale_fill_gradientn("Down-regulated\ngenes per MB", 
                         colours = down_col,
                         limits=c(0,max_count),
                         guide = guide_colorbar(order = 3)) +
    new_scale_fill() +
    geom_rect(data = up_df,
              aes(xmin=start, xmax=end, ymin=ymin, 
                  ymax=ymax, fill = gene_count)) +
    scale_fill_gradientn("Up-regulated\ngenes per MB",
                         colours = up_col,
                         limits=c(0,max_count),
                         guide = guide_colorbar(order = 1)) +
    geom_rect(data = chrom_df, aes(xmin =0, xmax = end, ymin = -1, ymax = 1),
              color = 'grey30', alpha = 0) +
    facet_grid(chrom~.) +
    scale_x_continuous(limits = c(0, max_chr_len * 10^deci_num), 
                       breaks = seq(0, max_chr_len * 10^deci_num, 10^deci_num),
                       labels = c(0, paste0(seq(1, max_chr_len), leading_zero, "M", sep = ""))) +
    ggtitle(plt_title) +
    theme(strip.text.y = element_text(size=rel(.8), 
                                      face="bold",
                                      angle = 0),
          axis.text.y = element_blank(),
          axis.ticks = element_blank(),
          panel.background = element_rect(fill = "white"),
          plot.title = element_text(hjust = 0.5))
  chrom_num <- length(chrom_order)
  p_height <- 2 * chrom_num / 7
  p_width <- max_chr_len * 8 / 3 
  ggsave(paste(out_prefix, 'png', sep='.'), 
         plot = p, width = p_width, height = p_height,
         dpi = 300, type = "cairo")
  ggsave(paste(out_prefix, 'pdf', sep='.'), 
         plot = p, width = p_width, height = p_height)
  
}


p <- arg_parser("for diff gene distribution")
p <- add_argument(p, "--diff_file", help = "diff analysis output")
p <- add_argument(p, "--chrom_size", help = "chr size")
p <- add_argument(p, "--out_dir", help = "output directory")
p <- add_argument(p, "--compare", help = "diff compare")
argv <- parse_args(p)

diff_file <- argv$diff_file
chrom_size <- argv$chrom_size
out_dir <- argv$out_dir
compare <- argv$compare

chrom_df <- read.delim(chrom_size, header = F, col.names = c('chrom', 'end'))
chrom_df$start <- 0
diff_df <- read.csv(diff_file)
diff_df$chrom <- factor(diff_df$chrom,
                        levels = chrom_df$chrom)
gene_type_num <- unique(diff_df$thickStart)
if (length(gene_type_num) > 1) {
  pcg_diff_df <- filter(diff_df, thickStart == 'protein_coding')
  pcg_prefix <- file.path(out_dir, paste('protein_coding', compare, 'diff_gene_location', sep = '_'))
  lnc_diff_df <- filter(diff_df, thickStart == 'lncRNA')
  lnc_prefix <- file.path(out_dir, paste('lncRNA', compare, 'diff_gene_location', sep = '_'))
  lnc_distribution_heatmap(pcg_diff_df, chrom_df, pcg_prefix, 1)
  lnc_distribution_heatmap(lnc_diff_df, chrom_df, lnc_prefix, 2)
} else {
  pcg_prefix <- file.path(out_dir, paste(compare, 'diff_gene_location', sep = '_'))
  lnc_distribution_heatmap(diff_df, chrom_df, pcg_prefix, 1)
}









