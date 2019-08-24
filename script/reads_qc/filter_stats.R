options(warn = -1)
options(stringsAsFactors = F)
suppressMessages(library('reshape2'))
suppressMessages(library('argparser'))
suppressMessages(library('plyr'))
suppressMessages(library('dplyr'))
suppressMessages(library('omplotr'))

p <- arg_parser('filter stats')
p <- add_argument(p,'--filter_dir',help = 'reads filter directory.')
argv <- parse_args(parser = p)


blank_theme <- theme_minimal()+
  theme(
    axis.title.x = element_blank(),
    axis.title.y = element_blank(),
    panel.border = element_blank(),
    panel.grid=element_blank(),
    axis.ticks = element_blank(),
    plot.title=element_text(size=14, face="bold"),
    legend.title = element_blank()
  )

reads_filter_plot <- function(filter_df, output) {
  sample_number <- length(unique(filter_df$sample))
  
  filter_df$percent <- paste(round(filter_df$proportion * 100, 3), '%', sep = '')
  filter_df$label <- paste(filter_df$stats,
                           '(',
                           filter_df$percent,
                           ')', sep = '')
  filter_df <- arrange(filter_df, desc(count))
  filter_df$stats <- factor(filter_df$stats,
                            levels = unique(filter_df$stats))
  stats_lab <- filter_df$label
  names(stats_lab) <- filter_df$stats
  stats_num <- length(filter_df$stats)
  stats_col <- RColorBrewer::brewer.pal(stats_num, 'Set2')
  names(stats_col) <- unique(filter_df$stats)
  plot_df <- filter_df[filter_df$count > 0, ]
  pie <- ggplot(plot_df, aes(x = "", y=proportion, fill = stats)) +
    geom_bar(width = 1, stat = "identity") +
    coord_polar("y", start = 0) +
    blank_theme + theme(axis.text.x=element_blank())
  
  if (sample_number > 1) {
    facet_wrap_ncol = round(sqrt(sample_number))
    pie <- pie + scale_fill_manual(values = stats_col) +
      facet_wrap(~sample, ncol = facet_wrap_ncol) +
      guides(fill=guide_legend(title = ''))
  } else {
    pie <- pie + scale_fill_manual(values = stats_col,
                      labels = stats_lab)
  }
  pie
  if (! is.null(output)) {
    plot_height <- 6 + sample_number/4
    plot_width <- 6 + sample_number/4
    save_ggplot(pie, output,
                width=plot_width,
                height=plot_height)
  }
  return(pie)
}

filter_dir <- argv$filter_dir

all_files <- list.files(filter_dir)
filter_files <- all_files[grep("*filter.txt", all_files)]


split_str <- function(strings, Split) {
  for (i in 1:nchar(strings)) {
    if (substr(strings, i, i) == Split) {
      return(c(substr(strings, 1, i - 1), substr(strings, i + 1, nchar(strings))))
    }
  }
}


file_list <- list()
for (i in seq(length(filter_files))) {
  each_sample_df <- read.delim(paste(filter_dir, filter_files[i], sep='/'),
                               header = F)
  colnames(each_sample_df) <- c('stats', 'count')
  each_sample_df$proportion <- prop.table(each_sample_df$count)
  sample_id <- split_str(filter_files[i], Split='.')[1]
  each_sample_df$sample <- sample_id
  each_sample_df[is.na(each_sample_df)] <- 0
  file_list[[i]] <- each_sample_df
  each_sample_out_name <- paste(sample_id, 'reads_filter', sep = '.')
  each_sample_out_path <- file.path(filter_dir, each_sample_out_name)
  reads_filter_plot(each_sample_df, each_sample_out_path)
}

#file_df <- ldply(file_list, data.frame)

#samples <- unique(file_df$sample)
#sample_number <- length(samples)
#selected_num <- ifelse(sample_number < 9, sample_number, 9)
#selected_df <- filter(file_df, sample %in% samples[1:selected_num])
#plot_out <- file.path(filter_dir, 'reads_filter.pie.report')
#selected_df$sample <- factor(selected_df$sample, levels = samples)
#reads_filter_plot(selected_df, plot_out)