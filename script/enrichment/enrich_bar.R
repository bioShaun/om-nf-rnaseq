suppressMessages(library(omplotr))
suppressMessages(library(stringr))
suppressMessages(library(argparser))

options(stringsAsFactors = F)
p <- arg_parser("enrichment bar plot")
p <- add_argument(p, "--enrich_file", help = "enrichment table file.")
argv <- parse_args(p)


wrap_long_name2 <- function(name, width=60) {
  if (nchar(name) > 120) {
      pre_name <- substr(name, 1, 60)
      last_name <- substr(name, nchar(name) - 54, nchar(name))
      name <- paste(pre_name, last_name, sep=' ... ')
  }
  return(paste(strwrap(name, width = width), collapse="\n"))
}


om_enrich_bar_plot2 <- function(enrich_df,
                               term_number=20,
                               ylab_title='-log10(qvalue)',
                               out_prefix=NULL,
                               title="") {
  enrich_df <- dplyr::filter(enrich_df, qvalue < 1)
  enrich_df <- head(enrich_df, term_number)
  if (dim(enrich_df)[1] > 0) {
    enrich_df$log10qvalue <- -log10(enrich_df$qvalue)
    enrich_df$color_val <- log(enrich_df$log10qvalue)
    enrich_df$wrap_term <- sapply(enrich_df$term, wrap_long_name2)
    enrich_df <- enrich_df[!duplicated(enrich_df$wrap_term), ]
    enrich_df$wrap_term <- factor(enrich_df$wrap_term,
                                  levels = rev(enrich_df$wrap_term))

    max_term_len <- max(nchar(as.character(enrich_df$term)))
    max_term_len <- ifelse(max_term_len > 60, 60, max_term_len)

    plot_witdh <- 6 + max_term_len *0.01
    plot_height <- dim(enrich_df)[1] / 3

    p <- ggplot(enrich_df, aes(wrap_term, log10qvalue, fill = ontology, alpha = color_val)) +
      geom_bar(stat = 'identity', width = 0.45, color='grey50') +
      coord_flip() +
      facet_grid(ontology~., scales = 'free_y',
                 space = 'free_y') +
      theme_onmath() +
      theme(axis.text.y = element_text(color = 'grey30', face = "bold",
                                       size = rel(0.75)),
            strip.text.y = element_text(angle = 0, colour = 'grey30')) +
      scale_fill_brewer(palette = 'Set1') +
      guides(fill = F, alpha = F) +
      xlab('') +
      ylab(ylab_title) + ggtitle(title)
    if (! is.null(out_prefix)) {
      save_ggplot(p, out_prefix,
                  width = plot_witdh,
                  height = plot_height)
    }
    return(p)
  } else {
    return('Nothing to plot!')
  }
}



enrich_file <- argv$enrich_file

out_name <- basename(enrich_file)
out_prefix <- str_replace(enrich_file, 'enrichment.txt', 'enrichment_barplot')
title <- str_replace(out_name, '_enrichment.txt', "")


enrich_plot_data_list <- clean_enrich_table(enrich_file)

om_enrich_bar_plot2(enrich_plot_data_list$table, 
                   ylab_title=enrich_plot_data_list$title, 
                   out_prefix=out_prefix, title=title)
