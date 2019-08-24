#!/usr/bin/env Rscript

# Command line argument processing
args <- commandArgs(trailingOnly=TRUE)

if (length(args) < 3) {
  stop("Usage: edgeR_heatmap_MDS.r <sample_1.bam> <sample_2.bam> <sample_3.bam> (more bam files optional)", call.=FALSE)
}

# Load / install required packages
if (!require("limma")){
    source("http://bioconductor.org/biocLite.R")
    biocLite("limma", suppressUpdates=TRUE)
    library("limma")
}
if (!require("edgeR")){
    source("http://bioconductor.org/biocLite.R")
    biocLite("edgeR", suppressUpdates=TRUE)
    library("edgeR")
}
if (!require("data.table")){
    install.packages("data.table", dependencies=TRUE, repos='http://cloud.r-project.org/')
    library("data.table")
}
if (!require("gplots")) {
    install.packages("gplots", dependencies=TRUE, repos='http://cloud.r-project.org/')
    library("gplots")
}

# Load count column from all files into a list of data frames
# Use data.tables fread as much much faster than read.table
# Row names are GeneIDs
temp <- lapply(lapply(args, fread, skip="Geneid", header=TRUE), function(x){return(as.data.frame(x)[,c(1, ncol(x))])})

# Merge into a single data frame
merge.all <- function(x, y) {
    merge(x, y, all=TRUE, by="Geneid")
}
data <- data.frame(Reduce(merge.all, temp))

# Clean sample name headers
colnames(data) <- gsub(".sorted.bam", "", colnames(data))

# Set GeneID as row name
rownames(data) <- data[,1]
data[,1] <- NULL

# Convert data frame to edgeR DGE object
dataDGE <- DGEList( counts=data.matrix(data) )

# Normalise counts
dataNorm <- calcNormFactors(dataDGE)

# Get the log counts per million values
logcpm <- cpm(dataNorm, prior.count=2, log=TRUE)
logcpm <- logcpm[rowSums(logcpm) > 0, ]
sample_cor <- cor(logcpm)
cor_dist <- as.dist(1-sample_cor)
hc = hclust(cor_dist)
order_sample_cor <- sample_cor[hc$order, hc$order]
write.csv(order_sample_cor, 'log2CPM_sample_cor.csv', quote=FALSE)