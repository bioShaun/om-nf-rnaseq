suppressMessages(library(dplyr))
suppressMessages(library(argparser))
suppressMessages(library(omplotr))

p <- arg_parser("perform go analysis")
p <- add_argument(p, "--name", help = "output name")
p <- add_argument(p, "--gene_list", help = "gene list to analysis")
p <- add_argument(p, "--go_anno", help = "gene id and transcript id mapping file")
p <- add_argument(p, "--gene_length", help = "gene length file")
p <- add_argument(p, "--out_dir", help = "diff analyssi output directory")
argv <- parse_args(p)



## get the parameter
gene_list <- argv$gene_list
name <- argv$name
go_anno_file <- argv$go_anno
gene_length_file <- argv$gene_length
go_out_dir <- argv$out_dir

## output directories
dir.create(go_out_dir, recursive = 1, showWarnings = F)
togo_dir <- file.path(go_out_dir, 'DAG')
dir.create(togo_dir, recursive = 1, showWarnings = F)

# read annotations
gene_length_df <- read.delim(gene_length_file, header = F)

# read gene list
gene_list_vector <- scan(gene_list, what=character())
# output file name
out_preifx_name <- paste(name, "go.enrichment", sep = ".")
out_prefix_path <- file.path(go_out_dir, out_preifx_name)
# run goseq
enrich_result <- om_goseq(gene_list_vector, gene_length_df,
  go_anno_file, out_prefix_path)
# run topgo
om_topgo(go_anno_file, gene_list_vector, enrich_result, name,
  togo_dir)
