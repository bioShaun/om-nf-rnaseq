suppressMessages(library(UpSetR))
suppressMessages(library(data.table))
suppressMessages(library(argparser))

options(stringsAsFactors = F)
p <- arg_parser("for upset plot")
p <- add_argument(p, "--gene_matrix", help = "a 1-0 matrix for plot")
p <- add_argument(p, "--out_prefix", help = "plot prefix")
argv <- parse_args(p)

gene_matrix <- argv$gene_matrix
out_prefix <- argv$out_prefix

exp_df <- fread(gene_matrix, check.names = F)
group_ids <- colnames(exp_df)[2:(ncol(exp_df))]

plot_height = 5 + length(group_ids) * 0.2
plot_width = 12 + max(nchar(group_ids)) / 20 

png(paste(out_prefix, "png", sep = "."),
        width = plot_width, height = plot_height,
        units = "in", res = 300)
upset(exp_df, sets = group_ids, sets.bar.color = "#56B4E9",
                        order.by = "freq", empty.intersections = "on",
                        nintersects = 50)
dev.off()

pdf(paste(out_prefix, "pdf", sep = "."),
        width = plot_width, height = plot_height, onefile = F)
upset(exp_df, sets = group_ids, sets.bar.color = "#56B4E9",
                        order.by = "freq", empty.intersections = "on",
                        nintersects = 50)
dev.off()
