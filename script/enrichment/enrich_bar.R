suppressMessages(library(omplotr))
suppressMessages(library(stringr))
suppressMessages(library(argparser))

options(stringsAsFactors = F)
p <- arg_parser("enrichment bar plot")
p <- add_argument(p, "--enrich_file", help = "enrichment table file.")
argv <- parse_args(p)


enrich_file <- argv$enrich_file

out_prefix <- str_replace(enrich_file, 'enrichment.txt', 'enrichment.barplot')

enrich_plot_data_list <- clean_enrich_table(enrich_file)

om_enrich_bar_plot(enrich_plot_data_list$table, 
                   ylab_title=enrich_plot_data_list$title, 
                   out_prefix=out_prefix)
