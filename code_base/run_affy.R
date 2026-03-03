#!/usr/bin/Rscript

# This script runs Affy using info read in from Comparison files with the following structure
# (tab-separated)
#
# ComparisonName      <outputs_prefix>
# Reference           <reference_conditions>
# file_path           condition
# <path_to_file_1>    <condition>
# ...                 ...
# <path_to_file_N>    <condition>

### ---------------------------------------- ###

parseArgs <- function() {
  
  # Read command line arguments
  args <- commandArgs()
  
  # Design file
  design_file <- args[match("--design", args) + 1]
  
  if("--norm_method" %in% args) { # RMA or GCRMA
    
    norm_method <- args[match("--norm_method", args) + 1]
    
  } else {
    
    norm_method <- "gcrma"
    
  }
  
  # ProbeID to Ensembl file
  probeid_to_ensembl <- args[match("--probeid_to_ensembl", args) + 1]
  
  if("--summarize_probes" %in% args) { # RMA or GCRMA
    
    summarize_probes <- TRUE
    
  } else {
    
    summarize_probes <- FALSE
    
  }
  
  return(c(design_file, norm_method, probeid_to_ensembl, summarize_probes))
  
}

### ---------------------------------------- ###

importDesign <- function(design_file) {
  
  # Reading design file
  temp <- read.delim(design_file, header = FALSE, sep = "\t", stringsAsFactors = FALSE)
  
  # Extracting analysis name
  analysis_name <- as.character(temp[1, 2])
  
  # Extracting reference
  reference <- as.character(temp[2, 2])
  
  # Generating sample_info dataframe
  sample_info <- as.data.frame(temp[4 : nrow(temp), 1 : ncol(temp)], stringsAsFactors = TRUE)
  colnames(sample_info) <- temp[3, 1 : ncol(temp)]
  
  return(list("analysis_name" = analysis_name, "reference" = reference, "sample_info" = sample_info))
  
}

### ---------------------------------------- ###

add_feature_data <- function(eset, file_path) {
  
  # Load probe_id to ensembl conversion table
  probeid_to_ensembl <- read.delim(file_path, sep="\t")
  probeid_to_ensembl <- subset(probeid_to_ensembl, probeid_to_ensembl$probe_id %in% rownames(expression_set))
  
  # Group table by probe_id, then keep probes mapping to only 1 gene
  probe_groups <- dplyr::group_by(probeid_to_ensembl, probe_id)
  probe_groups_summary <- dplyr::summarize(probe_groups, n_matches=n_distinct(ensembl_id))
  good_probes <- filter(probe_groups_summary, n_matches == 1)$probe_id
  
  # Filter eset
  eset <- eset[good_probes,]
  
  # Filter probeid_to_ensembl, sort as in rownames(eset), and add it to eset as feature_data
  probeid_to_ensembl <- subset(probeid_to_ensembl, probeid_to_ensembl$probe_id %in% good_probes)
  rownames(probeid_to_ensembl) <- probeid_to_ensembl$probe_id
  fData(eset) <- probeid_to_ensembl
  
  return(eset)
  
}

### ---------------------------------------- ###

merge_probes <- function(eset) {
  
  # Extract expression and feature data
  expression <- exprs(eset)
  feature_data <- fData(eset)
  
  # Group probes by genes
  gene_groups <- dplyr::group_by(feature_data, ensembl_id)
  gene_groups_summary <- dplyr::summarize(gene_groups, n_matches=n_distinct(probe_id))
  single_probe <- filter(gene_groups_summary, n_matches == 1)$ensembl_id
  multi_probe <- filter(gene_groups_summary, n_matches > 1)$ensembl_id
  
  # Average expression, then add genes with single probes
  averaged_expression <- lapply(multi_probe, FUN = function(gene) { colMeans(expression[subset(feature_data, feature_data$ensembl_id == gene)$probe_id,]) })
  averaged_expression <- as.data.frame(bind_rows(averaged_expression))
  rownames(averaged_expression) <- multi_probe
  
  # Extract expression of probes mapping to single genes
  unique_expression <- as.data.frame(expression[subset(feature_data, feature_data$ensembl_id %in% single_probe)$probe_id,])
  rownames(unique_expression) <- single_probe
  
  # Update expression set
  expression <- rbind(unique_expression, averaged_expression)
  feature_data <- data.frame(ensembl_id = rownames(expression),
                             row.names = rownames(expression))
  
  new_eset <- ExpressionSet(assayData = as.matrix(expression))
  fData(new_eset) <- feature_data
  
  return(new_eset)
  
}

### ---------------------------------------- ###

plot_pca <- function(expression, conditions, plot_name) {
  
  pca_data <- prcomp(expression, scale.=FALSE)
  explained_variance <- 100 * pca_data$sdev**2 / sum(pca_data$sdev**2)
  
  pca_plot_data <- data.frame(PC1 = pca_data$x[,1],
                              PC2 = pca_data$x[,2],
                              Condition = conditions)
  
  ggplot(pca_plot_data, aes(x=PC1, y=PC2, colour=Condition)) +
    geom_point(size=3) +
    xlab(paste("PC1 (", round(explained_variance[1], 2), ")%", sep="")) +
    ylab(paste("PC2 (", round(explained_variance[2], 2), ")%", sep="")) +
    theme(plot.title = element_text(hjust = 0.5)) +
    scale_shape_manual(values = c(4,15)) +
    scale_color_manual(values = c("darkorange2", "dodgerblue4"))
  ggsave(plot_name, dpi=300, width=6, height=5)
  
}

### ---------------------------------------- ###

plot_box <- function(expression, conditions, reference, plot_name) {
  
  boxcols <- rep("dodgerblue4", length(conditions))
  boxcols[conditions != reference] <- "darkorange2"
  #png(filename=plot_name, width=10, height=5)
  png(filename=plot_name)
  boxplot(expression, ann=FALSE, col=boxcols, axes=FALSE)
  dev.off()
  
}

### ---------------------------------------- ###

plot_volcano <- function(analysis, plot_name, p_thr=0.05) {
  
  # Volcano plot
  analysis <- subset(analysis, ! is.na(analysis$adj.P.Val))
  fold_change <- analysis$logFC
  log_pval <- - log10(analysis$adj.P.Val)
  genes_up <- sum(fold_change > 0 & log_pval > - log10(p_thr), na.rm = TRUE)
  genes_down <- sum(fold_change < 0 & log_pval > - log10(p_thr), na.rm = TRUE)
  genes_ns <- sum(log_pval <= - log10(p_thr), na.rm = TRUE)
  color <- rep(paste("NS (", genes_ns, ")", sep = ""), length(fold_change))
  color[fold_change < 0 & log_pval > - log10(as.double(p_thr))] <- paste("Down (", genes_down, ")", sep = "")
  color[fold_change > 0 & log_pval > - log10(as.double(p_thr))] <- paste("Up (", genes_up, ")", sep = "")
  color_palette <- c("green", "red", "gray")
  names(color_palette) <- c(paste("Up (", genes_up, ")", sep = ""),
                            paste("Down (", genes_down, ")", sep = ""),
                            paste("NS (", genes_ns, ")", sep = ""))
  qplot(fold_change, log_pval, fill = factor(color)) +
    geom_point(shape = 21, color = "black", size = 3, stroke = 0.5) +
    xlab("Log2 Fold Change") +
    ylab("- Log10 Pvalue") +
    scale_fill_manual(values = color_palette) +
    theme(legend.title = element_blank(),
          legend.text = element_text(size=14, face = "bold"),
          axis.text.y = element_text(size=12),
          axis.text.x = element_text(size=12),
          axis.title.y = element_text(size=14, face = "bold"),
          axis.title.x = element_text(size=14, face = "bold"),
          panel.background = element_blank(),
          panel.grid.major = element_blank(), 
          panel.grid.minor = element_blank(),
          panel.border = element_rect(colour = "black", fill=NA, linewidth=3))
  ggsave(filename=plot_name, width = 20, height = 25, units = "cm", dpi=300)
  
}

### ------------------MAIN------------------ ###

library(affy)
library(limma)
library(gcrma)
library(AnnotationDbi)
library(dplyr)
library(ggplot2)
library(RColorBrewer)
library(pheatmap)

parameters <- parseArgs()
print(parameters)

### Import design file

print("Reading design info")
design_info <- importDesign(parameters[1])

### Load CEL files

print("Loading CEL files")
cel_data <- ReadAffy(filenames=design_info$sample_info$file_path)

### Normalize data

if(parameters[2] == "gcrma") {
  
  print("GCRMA normalization")
  expression_set <- gcrma(cel_data)
  
} else {
  
  print("RMA normalization")
  expression_set <- rma(cel_data)
  
}

expression_set$sample <- colnames(expression_set)

### Add feature data

print("Adding feature data to expression set")
expression_set <- add_feature_data(expression_set, parameters[3])

### Average probes mapping to the same gene

if(parameters[4]) {
  
  print("Merging probes referring to the same gene")
  expression_set <- merge_probes(expression_set)
  
}

output_name <- paste(design_info$analysis_name, "_", parameters[2], "_counts.tsv", sep="")
write.table(exprs(expression_set), output_name, sep="\t", quote=FALSE)

### QC

# PCA
output_name <- paste(design_info$analysis_name, "_pca.png", sep="")
plot_pca(t(exprs(expression_set)), design_info$sample_info$condition, output_name)

# Boxplot
output_name <- paste(design_info$analysis_name, "_boxplot.png", sep="")
plot_box(exprs(expression_set), design_info$sample_info$condition, design_info$reference, output_name)

### Create design matrix

intercept <- rep(1, nrow(design_info$sample_info))
one_hot_conditions <- rep(0, nrow(design_info$sample_info))
one_hot_conditions[design_info$sample_info$condition != design_info$reference] <- 1

design_matrix <- as.matrix(data.frame(intercept = intercept,
                                      contrast = one_hot_conditions,
                                      row.names = basename(design_info$sample_info$file_path)))

### Fit linear model

# Fit model
print("Fitting model")
model_fit <- lmFit(expression_set, design_matrix)
#model_fit <- eBayes(model_fit, robust=TRUE)
model_fit <- eBayes(model_fit)

# Extract results
model_fit_results <- limma::topTable(model_fit, coef="contrast", adjust="BH", p.value=1, lfc=0, number=nrow(model_fit))

# Write to file
output_name <- paste(design_info$analysis_name, "_affy.tsv", sep="")
write.table(model_fit_results, output_name, sep="\t", quote=FALSE, row.names=FALSE)

# Volcano plot
output_name <- paste(design_info$analysis_name, "_volcano.png", sep="")
plot_volcano(model_fit_results, output_name, 0.05)
