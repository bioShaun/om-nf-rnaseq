suppressMessages(library(ggplot2))
suppressMessages(library(dplyr))
suppressMessages(library(RColorBrewer))
suppressMessages(library(scales))
suppressMessages(library(argparser))
suppressMessages(library(stringr))

options(stringsAsFactors = F)

p <- arg_parser("for diff gene distribution")
p <- add_argument(p, "--lnc_cis", help = "lnc neighbor gene file")
p <- add_argument(p, "--chrom_size", help = "chr size")
p <- add_argument(p, "--out_dir", help = "output directory")
argv <- parse_args(p)


lnc_cis <- argv$lnc_cis
chrom_size <- argv$chrom_size
out_dir <- argv$out_dir

lnc_cis_df <- read.csv(lnc_cis)
lnc_cis_df$pos <- (lnc_cis_df$start + lnc_cis_df$end) / 2

chr_df <- read.delim(chrom_size, header = F, col.names = c('chrom', 'end'))
chr_order <- chr_df$chrom


lnc_cis_df <- filter(lnc_cis_df, chrom %in% chr_order)
chr1_lnc_cis_df <- filter(lnc_cis_df, chrom == '1')


pos_pcc_col <- colorRampPalette(brewer.pal(9, 'Reds'))(100)
neg_pcc_col <- rev(colorRampPalette(brewer.pal(9, 'Blues'))(100))
pcc_col <- c(neg_pcc_col, pos_pcc_col)
max_pcc <- max(abs(lnc_cis_df$pcc))

lnc_cis_df$log_count <- log2(lnc_cis_df$gene_count + 1)

deci_num <- floor(log10(max(chr_df$end)))
max_chr_len = ceiling(max(chr_df$end)/(10^deci_num))
leading_zero <- paste(rep("0", deci_num - 6), collapse = "")
chr_df$chrom <- factor(chr_df$chrom,
                         levels = chr_order)
lnc_cis_df$chrom <- factor(lnc_cis_df$chrom,
                       levels = chr_order)
p <- ggplot() +
  geom_segment(data = chr_df, aes(x=0, xend=end, y=0, yend=0), color = 'grey70') +
  geom_rect(data = lnc_cis_df,
            aes(xmin=start, xmax=end, ymin=0, 
                ymax=log_count, fill = pcc)) +
  scale_fill_gradientn("PCC", 
                       colours = pcc_col,
                       limits=c(-max_pcc,max_pcc),
                       guide = guide_colorbar(order = 3)) +
  facet_grid(chrom~.) +
  theme(strip.text.y = element_text(size=rel(.8),
                                    face="bold",
                                    angle = 0),
        axis.text.y = element_blank(),
        axis.ticks = element_blank(),
        panel.background = element_rect(fill = "white"),
        plot.title = element_text(hjust = 0.5)) +
  ylab('') + xlab('') + ggtitle('lncRNA PCG Neighbour Distribution') +
  scale_x_continuous(limits = c(0, max_chr_len * 10^deci_num), 
                     breaks = seq(0, max_chr_len * 10^deci_num, 10^deci_num),
                     labels = c(0, paste0(seq(1, max_chr_len), leading_zero, "M", sep = "")))
  
chrom_num <- length(chr_order)
p_height <- 2 * chrom_num / 7
p_width <- max_chr_len * 8 / 3 

out_prefix <- file.path(out_dir, 'lncRNA_PCG_neighbour_distribution')
ggsave(paste(out_prefix, 'png', sep='.'), 
       plot = p, width = p_width, height = p_height,
       dpi = 300, type = "cairo")
ggsave(paste(out_prefix, 'pdf', sep='.'), 
       plot = p, width = p_width, height = p_height)

