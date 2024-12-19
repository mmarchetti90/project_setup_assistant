#!/usr/bin/Rscript

# This script runs DEXSeq using information read in from Comparison files
# N.B. This is designed for simple binary comparisons as described in a comparison design file in the column "condition" and field "reference"

### ---------------------------------------- ###

parseArgs <- function() {

	# Read command line arguments
	args <- commandArgs()

	# Counts file
	counts_file <- args[match("--counts", args) + 1]
	
	# Design file
	design_file <- args[match("--design", args) + 1]
	
	# Threads
	threads <- args[match("--threads", args) + 1]

	return(c(counts_file, design_file, threads))

}

### ---------------------------------------- ###

importDesign <- function(design_file) {

	# Reading design file
	temp <- read.delim(design_file, header = FALSE, sep = "\t", stringsAsFactors = FALSE)

	# Extracting analysis name
	analysis_name <- as.character(temp[1, 2])

	# Extracting comparison formula
	if(! "~ sample + exon + " %in% temp[2, 2]) {
	  
	  condition <- as.formula(paste("~ sample + exon + ", temp[4, 2], ":exon", sep=""))
	  
	} else {
	  
	  condition <- as.formula(temp[2, 2])
	  
	}
	
	# Extracting reference
	reference <- as.character(temp[3, 2])

	# Generating sample_info dataframe
	sample_info <- as.data.frame(temp[5 : nrow(temp), 2 : ncol(temp)], stringsAsFactors = TRUE)
	rownames(sample_info) <- temp[5 : nrow(temp), 1]
	colnames(sample_info) <- temp[4, 2 : ncol(temp)]

	return(list("analysis_name" = analysis_name, "formula" = condition, "reference" = reference, "sample_info" = sample_info))

}

### ---------------------------------------- ###

runDEA <- function(params, cnts, des) {

	# Removing samples from design not in counts table
  	des$sample_info <- subset(des$sample_info, rownames(des$sample_info) %in% colnames(cnts))
  
  	# Setting number of threads
  	threads <- MulticoreParam(as.integer(params[3]))
  
  	# Creating a DEXSeq data matrix
  	dxd <- DEXSeqDataSet(countData = cnts[, rownames(des$sample_info)], sampleData = des$sample_info, design = des$formula, featureID = cnts$ExonID, groupID = cnts$GeneID)
  
  	# Normalization
  	dxd <- estimateSizeFactors(dxd)
  
  	# Estimate data dispersion
  	dxd <- estimateDispersions(dxd, BPPARAM = threads, quiet = FALSE)
  
  	# Differential Exon Usage (DEU) testing
  	dxd <- testForDEU(dxd, BPPARAM = threads)
  
  	dxd <- estimateExonFoldChanges(dxd, BPPARAM = threads, fitExpToVar = "condition", denominator = des$reference)
  
  	# Save RDS object
  	saveRDS(dxd, file = paste("DEU_", des$analysis_name, ".rds", sep = ""))
  
  	# Extracting results
  	analysis <- DEXSeqResults(dxd)
  	analysis_dataframe <- as.data.frame(analysis)
  	analysis_dataframe <- analysis_dataframe[order(analysis_dataframe$padj, decreasing = FALSE),
                                           colnames(analysis_dataframe) != "genomicData"]
  	output_name <- paste("DEU_", des$analysis_name, ".tsv", sep = "")
  	write.table(analysis_dataframe, output_name, row.names = FALSE, sep = "\t")
  
  	### QUALITY CONTROL PLOTS ---------------- ###
  
  	# Plot dispersion estimate
  	output_name <- paste("DEU_", des$analysis_name, "_Dispersion-Plot.png", sep = "")
  	png(file = output_name, width = 600, height = 600)
  	plotDispEsts(dxd)
  	dev.off()
  
  	# MA plot
  	output_name <- paste("DEU_", des$analysis_name, "_MA-Plot.png", sep = "")
  	png(file = output_name, width = 600, height = 600)
  	plotMA(analysis, ylim=c(-5, 5))
  	dev.off()
  
}

### ------------------MAIN------------------ ###

library(DEXSeq)
library(BiocParallel)

parameters <- parseArgs()

# Import counts and experimental design
counts <- read.delim(as.character(parameters[1]), header = TRUE, sep = "\t", check.names = FALSE)
design <- importDesign(as.character(parameters[2]))

# Running DESeq2
runDEA(parameters, counts, design)
