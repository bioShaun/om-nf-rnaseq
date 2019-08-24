library(pathview)

args<-commandArgs(T)
spe <- args[1]
pathway_id <- args[2]
spe_gene <- args[3]
database <- args[4]
out_dir <- args[5]

options(bitmapType='cairo')

database_dir <- paste(database, spe, sep = '/')

pathway_data <- read.table(spe_gene, header = F, sep = '\t', row.names = 1)

setwd(out_dir)

pv.out <- pathview(gene.data = pathway_data, gene.idtype="kegg",
                   pathway.id = pathway_id, species = spe,
                   kegg.dir = database_dir,min.nnodes = 1,
                   kegg.native = T, same.layer=T,
                   low = list(gene = "SteelBlue", cpd = "blue"),
                   mid = list(gene = "LemonChiffon", cpd = "gray"),
                   high = list(gene = "OrangeRed", cpd = "yellow"))

pv.out <- pathview(gene.data = pathway_data, gene.idtype="kegg",
                   pathway.id = pathway_id, species = spe,
                   kegg.dir = database_dir, min.nnodes = 1 ,
                   out.suffix = "pathview.gene",
                   kegg.native = T, same.layer=F,
                   low = list(gene = "SteelBlue", cpd = "blue"),
                   mid = list(gene = "LemonChiffon", cpd = "gray"),
                   high = list(gene = "OrangeRed", cpd = "yellow"))

