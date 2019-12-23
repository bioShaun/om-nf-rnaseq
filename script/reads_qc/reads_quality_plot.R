#2017-03-01
options(warn = -1)
options(stringsAsFactors = F)
suppressMessages(library('reshape2'))
suppressMessages(library('argparser'))
suppressMessages(library('plyr'))
suppressMessages(library('dplyr'))
suppressMessages(library('omplotr'))

p <- arg_parser('reads quality plot')
p <- add_argument(p,'--rq_dir',help = 'reads quality stats file directory.')
argv <- parse_args(parser = p)

# read parameters
rq_dir <- argv$rq_dir

all_files <- list.files(rq_dir)
rq_files <- all_files[grep("*reads_quality.txt", all_files)]


split_str <- function(strings, Split) {
  for (i in 1:nchar(strings)) {
    if (substr(strings, i, i) == Split) {
      return(c(substr(strings, 1, i - 1), substr(strings, i + 1, nchar(strings))))
    }
  }
}


reads_quality_plot2 <- function(plot_data, output=NULL, title="") {
  
  sample_number <- length(unique(plot_data$sample))
  
  max_qual <- max(plot_data$quality)
  p <- ggplot(plot_data, aes(X.Base, quality)) +
    geom_bar(width = 0.75, stat = 'identity',
             fill = 'lightblue') +
    theme_onmath() +
    scale_x_continuous(breaks = seq(0, 300, 50)) +
    geom_vline(xintercept = 150, 
               color = 'grey50',
               lty = 2, size=0.5) +
    xlab("Postion") + ylab("Quality") +
      coord_cartesian(ylim=c(20,max_qual)) +
      ggtitle(title)
  
  if (sample_number > 1) {
    facet_wrap_ncol = round(sqrt(sample_number))
    p <- p + facet_wrap(~sample, ncol = facet_wrap_ncol)
  }
  
  
  
  if (! is.null(output)) {
    plot_height <- 6 + sample_number/4
    plot_width <- 8 + sample_number/4
    save_ggplot(p, output,
                width=plot_width,
                height=plot_height)
  }
  return(p)
}

rq_file_list <- list()
for (i in seq(length(rq_files))) {
  each_sample_rq_df <- read.delim(paste(rq_dir, rq_files[i], sep='/'))
  sample_id <- split_str(rq_files[i], Split='.')[1]
  each_sample_rq_df$sample <- sample_id
  each_sample_rq_df[is.na(each_sample_rq_df)] <- 0
  rq_file_list[[i]] <- each_sample_rq_df
  each_sample_out_name <- paste(sample_id, 'reads_quality.bar', sep = '.')
  each_sample_out_path <- file.path(rq_dir, each_sample_out_name)
  reads_quality_plot2(each_sample_rq_df, each_sample_out_path, title=sample_id)
}

#rq_file_df <- ldply(rq_file_list, data.frame)

#samples <- unique(rq_file_df$sample)
#sample_number <- length(samples)
#selected_num <- ifelse(sample_number < 9, sample_number, 9)
#selected_df <- filter(rq_file_df, sample %in% samples[1:selected_num])
#rq_plot_out <- file.path(rq_dir, 'reads_quality.bar.report')
#selected_df$sample <- factor(selected_df$sample, levels = samples)
#reads_quality_plot2(selected_df, rq_plot_out)
