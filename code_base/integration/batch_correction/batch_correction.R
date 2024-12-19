#!/usr/bin/Rscript

# This script runs ComBat_seq on a counts matrix

# INPUT
# -----
# counts
#   Counts matrix with shape M * (N + 1), with M being features and N samples.
#   The first column should be features unique identifiers
# 
# sample_info
#   Table with two required columns:
#     Sample -> unique sample identifiers (must match count matrix columns)
#     Batch -> integer defining the batch of each sample

### ---------------------------------------- ###

parseArgs <- function() {

	# Read command line arguments
	args <- commandArgs()

	# Counts file
	counts_file <- args[match("--counts", args) + 1]
	
	# Sample info file
	sample_info_file <- args[match("--sample_info", args) + 1]

	return(c(counts_file, sample_info_file))

}

### ---------------------------------------- ###

batchCorrection <- function(file_name, cnts, info) {
  
  # Removing samples not in the matrix
  info <- subset(info, info$Sample %in% colnames(cnts))
  
  # Correcting for batch effect using ComBat_seq
  cnts[, info$Sample] <- ComBat_seq(as.matrix(cnts[, info$Sample]),
                                    info$Batch,
                                    group = NULL,
                                    full_mod = TRUE)
  cnts <- as.data.frame(cnts)
  
  # Save to file
  out_name <- paste(substr(file_name, 1, nchar(file_name) - 4), "_BatchCorrected.tsv", sep="")
  write.table(cnts, out_name, quote=FALSE, sep='\t', row.names=FALSE, col.names=TRUE)
  
}

### ------------------MAIN------------------ ###

library(sva)

parameters <- parseArgs()

# Import counts
counts <- read.delim(as.character(parameters[1]), header = TRUE, sep = "\t", check.names = FALSE)

# Import sample info
sample_info <- read.delim(as.character(parameters[2]), header = TRUE, sep = "\t", check.names = FALSE)

# Batch correction
batchCorrection(as.character(parameters[1]), counts, sample_info)
