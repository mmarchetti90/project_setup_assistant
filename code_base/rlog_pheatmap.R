#!/usr/bin/Rscript

# Plots RLOG heatmaps downstream of run_deseq2.R

### ---------------------------------------- ###

parseArgs <- function() {
  
  # Read command line arguments
  args <- commandArgs()
  
  # Load RLOG file
  rlog_file <- args[match("--rlog_file", args) + 1]
  
  # Load DEA file(s)
  dea_files <- args[match("--dea_files", args) + 1]
  
  # Gene list file location
  gene_list_file <- args[match("--gene_list_file", args) + 1]
  
  # Gene list name
  gene_list_name <- args[match("--gene_list_name", args) + 1]

  # Samples file
  samples_file <- args[match("--samples_file", args) + 1]
  
  return(c(rlog_file, dea_files, gene_list_file, gene_list_name, samples_file))
  
}

### ---------------------------------------- ###

plotHeatmap <- function(out, data, ann_row, ann_col, ann_colors, c_rows=TRUE, c_cols=FALSE) {
  
  heatmap_colors <- colorRampPalette(c("blue", "white", "red"))(50)

  res <- pheatmap(data,
                  border_color='black',
                  color=heatmap_colors,
                  cellwidth=30,
                  cellheight=15,
                  scale="row",
                  cluster_rows=c_rows,
                  cluster_cols=c_cols,
                  annotation_row = ann_row,
                  annotation_col = ann_col,
                  annotation_colors = ann_colors,
                  show_colnames = FALSE,
                  legend=TRUE,
                  angle_col=90,
                  fontsize=12,
                  breaks=seq(from=-2, to=2, length.out=51),
                  filename=out)

  return(res)
  
}

### ------------------MAIN------------------ ###

library(pheatmap)
library(RColorBrewer)

# Read args
parameters <- parseArgs()
print(parameters)

# Load data
rlog <- read.delim(parameters[1], header = TRUE, row.names = 1, check.names = FALSE, sep = "\t")
genes_of_interest <- read.delim(parameters[3], row.names = 1, sep="\t", check.names = FALSE)
samples_table <- read.delim(parameters[5], row.names = 1, sep="\t", check.names = FALSE)

p_thr <- 0.05
dea_files <- unlist(strsplit(parameters[2], ','))
dea_data <- list()
for(df in dea_files) {
  
  analysis_name <- gsub("_annotated.tsv", "", gsub("DEA_", "", basename(df)))
  data <- read.delim(df, sep = "\t", check.names = FALSE)
  data <- subset(data, data$padj < p_thr)
  dea_data[[analysis_name]] <- data
  
}

# Filter for genes of interest
genes_of_interest <- subset(genes_of_interest, rownames(genes_of_interest) %in% rownames(rlog))
rlog_sub <- rlog[rownames(genes_of_interest),]
rownames(rlog_sub) <- genes_of_interest$GeneSymbol

# Order samples by type
rlog_sub <- rlog_sub[, rownames(samples_table)]

# Row annotation matrix
genes <- genes_of_interest$GeneSymbol
ann_row <- list()
for(name in names(dea_data)) {
  
  upreg <- subset(dea_data[[name]], dea_data[[name]]$log2FoldChange > 0)$gene_symbol
  downreg <- subset(dea_data[[name]], dea_data[[name]]$log2FoldChange < 0)$gene_symbol
  ann <- rep('Ns', length(genes))
  ann[genes %in% upreg] <- "Upreg"
  ann[genes %in% downreg] <- "Downreg"
  ann <- factor(ann, levels = c("Downreg", "Ns", "Upreg"))
  ann_row[[name]] <- ann
  
}
ann_row <- data.frame(ann_row, row.names = genes)

# Init annotation colors
ann_colors <- list()

# Add conditions colors
ann_colors[["condition"]] <- c()
conditions <- unique(samples_table$condition)
for(i in 1:length(conditions)) {
  
  color = c("#ff9289", "#96ca00", "#00dae0", "#e199ff")[i] # Max 4 conditions in datasets
  ann_colors[["condition"]] <- c(ann_colors[["condition"]], color)
  
}

names(ann_colors[["condition"]]) <- conditions

# Add row colors
for(col in colnames(ann_row)) {
  
  ann_colors[[col]] <- c()
  
  if("Downreg" %in% ann_row[[col]]) {
    
    ann_colors[[col]] <- c(ann_colors[[col]], c("Downreg" = "blue"))
    
  }
  
  if("Ns" %in% ann_row[[col]]) {
    
    ann_colors[[col]] <- c(ann_colors[[col]], c("Ns" = "white"))
    
  }
  
  if("Upreg" %in% ann_row[[col]]) {
    
    ann_colors[[col]] <- c(ann_colors[[col]], c("Upreg" = "red"))
    
  }
  
}

# Reverse order of annotations to match sample order
if(ncol(ann_row) > 1) {
  
  ann_row <- ann_row[, rev(colnames(ann_row))]
  
}

# Heatmaps
output_prefix <- parameters[4]
output_name <- paste(output_prefix, "_rlog-heatmap_hclust_allgenes.pdf", sep='')
res <- plotHeatmap(output_name, rlog_sub, ann_row, samples_table, ann_colors, TRUE, FALSE)

# Order genes according to clustering
#row_order <- res$tree_row[["order"]]
#rlog_sub <- rlog_sub[row_order,]

# Split Heatmaps, keeping same order as full one
#output_name <- paste(output_prefix, "_rlog-heatmap_sub1.pdf", sep='')
#null_res <- plotHeatmap(output_name, rlog_sub[1:floor(nrow(rlog_sub) / 2),], ann_row, samples_table, ann_colors, FALSE, FALSE)

#output_name <- paste(output_prefix, "_rlog-heatmap_sub2.pdf", sep='')
#null_res <- plotHeatmap(output_name, rlog_sub[(floor(nrow(rlog_sub) / 2) + 1) : nrow(rlog_sub),], ann_row, samples_table, ann_colors, FALSE, FALSE)

# Remove genes that are not significant in any condition
good_genes <- names(which(rowSums(ann_row == "Ns") < ncol(ann_row)))
ann_row <- subset(ann_row, rownames(ann_row) %in% good_genes)
rlog_sub <- rlog_sub[good_genes,]

# Heatmaps
output_prefix <- parameters[4]
output_name <- paste(output_prefix, "_rlog-heatmap_hclust_significantgenes.pdf", sep='')
res <- plotHeatmap(output_name, rlog_sub, ann_row, samples_table, ann_colors, TRUE, FALSE)
