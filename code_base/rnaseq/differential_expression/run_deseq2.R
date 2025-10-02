#!/usr/bin/Rscript

# This script runs DESeq2 using information read in from Comparison files
# N.B. Can run both Wald or LRT

### ---------------------------------------- ###

parseArgs <- function() {

	# Read command line arguments
	args <- commandArgs()

	# Counts file
	counts_file <- args[match("--counts", args) + 1]
	
	# Design file
	design_file <- args[match("--design", args) + 1]
	
	# Gene counts lower limit
	min_reads <- args[match("--mincounts", args) + 1]

	# Adjusted p value threshold
	pval <- args[match("--p_thr", args) + 1]
	
	# DESeq2 mode, i.e. Wald test or LRT
	if("--mode" %in% args) {
	  
	  run_mode <- args[match("--mode", args) + 1]
	  
	} else {
	  
	  run_mode <- "wald"
	  
	}

	return(c(counts_file, design_file, min_reads, pval, run_mode))

}

### ---------------------------------------- ###

importDesign <- function(design_file, design_type) {

	# Reading design file
	temp <- read.delim(design_file, header = FALSE, sep = "\t", stringsAsFactors = FALSE)
	
	if(design_type == "wald") {
	  
	  # Extracting analysis name
	  analysis_name <- as.character(temp[1, 2])
	  
	  # Extracting comparison formula
	  condition <- as.formula(temp[2, 2])
	  
	  # Extracting reference
	  reference <- as.character(temp[3, 2])
	  
	  # Generating sample_info dataframe
	  sample_info <- as.data.frame(temp[5 : nrow(temp), 2 : ncol(temp)], stringsAsFactors = TRUE)
	  rownames(sample_info) <- temp[5 : nrow(temp), 1]
	  colnames(sample_info) <- temp[4, 2 : ncol(temp)]
	  
	  return(list("analysis_name" = analysis_name, "full_model" = condition, "reference" = reference, "sample_info" = sample_info))
	  
	} else {
	  
	  # Extracting analysis name
	  analysis_name <- as.character(temp[1, 2])
	  
	  # Extracting full model formula
	  full_model <- as.formula(temp[2, 2])
	  
	  # Extracting reduced model formula
	  reduced_model <- as.formula(temp[3, 2])
	  
	  # Extracting reference
	  reference <- as.character(temp[4, 2])
	  
	  # Generating sample_info dataframe
	  sample_info <- as.data.frame(temp[6 : nrow(temp), 2 : ncol(temp)], stringsAsFactors = TRUE)
	  rownames(sample_info) <- temp[6 : nrow(temp), 1]
	  colnames(sample_info) <- temp[5, 2 : ncol(temp)]
	  
	  return(list("analysis_name" = analysis_name, "full_model" = full_model, "reduced_model" = reduced_model, "reference" = reference, "sample_info" = sample_info))
	  
	}

}

### ---------------------------------------- ###

runDEA <- function(params, cnts, des) {

	# Removing samples from design not in counts table
  des$sample_info <- subset(des$sample_info, rownames(des$sample_info) %in% colnames(cnts))
  	
	# Filtering and ordering counts samples
	cnts <- cnts[, rownames(des$sample_info)]

	# Creating a DESeq2 data matrix
	dds <- DESeqDataSetFromMatrix(countData = cnts, colData = des$sample_info, design = des$full_model)

	# Filtering out genes with too few reads in more than half the samples
	dds <- dds[rowSums(counts(dds) >= as.integer(params[3])) >= ncol(dds) / 2,]

	# Finding parameter to relevel
	for(col in colnames(des$sample_info)) {
	  
		if(des$reference %in% des$sample_info[,col]) {
		  
			comparison_parameter <- col
			break
			
		}
	  
	}

	# Releveling parameter to reference
	dds[[comparison_parameter]] <- relevel(dds[[comparison_parameter]], ref=des$reference)
	dds[[comparison_parameter]] <- droplevels(dds[[comparison_parameter]])

	# Differential expression analysis
	if(params[5] == "wald") {
	  
	  dds <- DESeq(dds, test = "Wald", parallel = T)
	  
	  contrast <- c(comparison_parameter, rev(levels(dds[[comparison_parameter]])))
	  analysis <- results(dds, alpha = as.numeric(params[4]), contrast = contrast)
	  #analysis <- lfcShrink(dds, contrast = contrast, res = analysis, type = 'ashr')
	  
	} else {
	  
	  dds <- DESeq(dds, test = "LRT", reduced = des$reduced_model, parallel = T)
	  
	  analysis <- results(dds, alpha = as.numeric(params[4]))
	  #analysis <- lfcShrink(dds, contrast = contrast, res = analysis, type = 'ashr')
	  
	}
	
	# Save RDS object and results table
	output_name <- paste("DEA_", des$analysis_name, ".rds", sep = "")
	saveRDS(dds, file = output_name)

	output_name <- paste("DEA_", des$analysis_name, ".tsv", sep = "")
	write.table(analysis, output_name, sep = "\t")

	### QUALITY CONTROL PLOTS ---------------- ###
	
	# MA plot
	output_name <- paste("DEA_", des$analysis_name, "_MA-Plot.png", sep = "")
	png(file = output_name, width = 600, height = 600)
	plotMA(analysis, ylim=c(-5, 5))
	dev.off()

	# Regularized log transformation
	rld <- rlog(dds, blind = T)
	
	# Distance heatmap
	sampleDists <- dist(t(assay(rld)))
	colors <- colorRampPalette(rev(brewer.pal(n = 9, name = "Blues")))(255)
	output_name <- paste("DEA_", des$analysis_name, "_SampleDistance.png", sep = "")
	pheatmap(as.matrix(sampleDists), clustering_distance_rows = sampleDists, clustering_distance_cols = sampleDists, color = colors, filename = output_name)
	
	# PCA
	output_name <- paste("DEA_", des$analysis_name, "_PCA.png", sep = "")
	png(file = output_name, width = 600, height = 600)
	pca_plot <- plotPCA(rld, intgroup = c(comparison_parameter))
	print(pca_plot)
	dev.off()
	
	# Plotting variance vs mean value
	gene_mean <- log10(rowMeans(assay(rld)))
	gene_variance <- log10(rowVars(assay(rld)))
	output_name <- paste("DEA_", des$analysis_name, "_MeanVsVariance.png", sep = "")
	qplot(gene_mean, gene_variance) +
	  geom_point() +
	  xlab("Log10 Gene Mean") +
	  ylab("Log10 Gene Variance")
	ggsave(filename=output_name, dpi=300)

	# Count genes by behaviour
	genes_up <- sum(analysis$log2FoldChange > 0 & analysis$padj < as.double(params[4]), na.rm = TRUE)
	genes_down <- sum(analysis$log2FoldChange < 0 & analysis$padj < as.double(params[4]), na.rm = TRUE)
	genes_ns <- sum(analysis$padj > as.double(params[4]), na.rm = TRUE)

	# Remove genes with NA padj
	analysis <- subset(analysis, ! is.na(analysis$padj))

	# Prepare data for Volcano plot
	fold_change <- analysis$log2FoldChange
	log_pval <- - log10(analysis$padj)
	color <- rep(paste("NS (", genes_ns, ")", sep = ""), length(fold_change))
	color[fold_change < 0 & log_pval > - log10(as.double(params[4]))] <- paste("Down (", genes_down, ")", sep = "")
	color[fold_change > 0 & log_pval > - log10(as.double(params[4]))] <- paste("Up (", genes_up, ")", sep = "")
	color_palette <- c("green", "red", "gray")
	names(color_palette) <- c(paste("Up (", genes_up, ")", sep = ""),
	                          paste("Down (", genes_down, ")", sep = ""),
	                          paste("NS (", genes_ns, ")", sep = ""))
	
	# Volcano plot
	output_name <- paste("DEA_", des$analysis_name, "_VolcanoPlot.png", sep = "")
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
	        panel.border = element_rect(colour = "black", fill=NA, size=3))
	ggsave(filename=output_name, width = 20, height = 25, units = "cm", dpi=300)

}

### ------------------MAIN------------------ ###

library(DESeq2)
library(ggplot2)
library(RColorBrewer)
library(pheatmap)

# Read command line parameters
parameters <- parseArgs()

# Import counts and experimental design
counts <- read.delim(as.character(parameters[1]), header = TRUE, row.names = 1, sep = "\t", check.names = FALSE)
design <- importDesign(as.character(parameters[2]), parameters[5])

# Running DESeq2
runDEA(parameters, counts, design)
