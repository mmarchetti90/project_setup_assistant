#!/usr/bin/Rscript

# This script converts a raw counts matrix to RLOG using DESeq2
# Input counts table should be tab-separated with the first column reporting unique gene
# identifiers and the following should report the gene counts for individual samples

### ---------------------------------------- ###

parseArgs <- function() {

	# Read command line arguments
	args <- commandArgs()

	# Counts file
	counts_file <- args[match("--counts", args) + 1]
	
	# Design file
	design_file <- args[match("--design", args) + 1]

	return(c(counts_file, design_file))

}

### ---------------------------------------- ###

importDesign <- function(design_file) {

	# Reading design file
	temp <- read.delim(design_file, header = FALSE, sep = "\t", stringsAsFactors = FALSE)

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

	return(list("analysis_name" = analysis_name, "formula" = condition, "reference" = reference, "sample_info" = sample_info))

}

### ---------------------------------------- ###

getRlog <- function(cnts, des) {
  
	# Removing samples from design not in counts table
	des$sample_info <- subset(des$sample_info, rownames(des$sample_info) %in% colnames(cnts))
  	
	# Filtering and ordering counts samples
	cnts <- cnts[, rownames(des$sample_info)]

	# Creating a DESeq2 data matrix
	dds <- DESeqDataSetFromMatrix(countData = cnts, colData = des$sample_info, design = des$formula)

	# Regularized log transformation
	#rld <- rlog(dds, blind = T)
	rld <- vst(dds, blind = T)

	# Convert to data frame
	rld <- assay(rld)

	rld <- cbind(data.frame('GeneID' = rownames(rld)),
          rld)
	
	return(rld)

}

### ------------------MAIN------------------ ###

library(DESeq2)

parameters <- parseArgs()

# Import counts and experimental design
counts <- read.delim(as.character(parameters[1]), header = TRUE, row.names = 1, sep = "\t", check.names = FALSE)
design <- importDesign(as.character(parameters[2]))

# Rlog transformation
rlog_counts <- getRlog(counts, design)

# Save to file
write.table(rlog_counts, "rlog_counts.tsv", sep="\t", quote=FALSE, row.names=FALSE)
