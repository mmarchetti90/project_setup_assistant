#!/usr/bin/Rscript

# Uses the MergedGeneCounts.tsv file generate by gdc_merge_counts.py and generates DESeq2 RLOG counts

### ---------------------------------------- ###

getRlog <- function(cnts) {
  
  # Create smple info table as colData
  sample_info <- data.frame(sample = rep(1, ncol(cnts) - 2),
                            row.names = colnames(cnts)[3 : ncol(cnts)])
  
  # Record gene names as row.names
  rownames(cnts) <- cnts$GeneID
  gene_info <- cnts[, 1:2]
  cnts <- cnts[3:ncol(cnts)]

	# Creating a DESeq2 data matrix
	dds <- DESeqDataSetFromMatrix(countData = cnts, colData = sample_info, design = ~ 1)

	# Regularized log transformation
	#rld <- rlog(dds, blind = T)
	rld <- vst(dds, blind = T)

	# Convert to data frame
	rld <- assay(rld)

	# Add gene info
	rld <- cbind(gene_info, rld)
	
	return(rld)

}

### ------------------MAIN------------------ ###

library(DESeq2)

# Import counts and experimental design
counts <- read.delim("MergedGeneCounts.tsv", header = TRUE, row.names = NULL, sep = "\t", check.names = FALSE)

# Rlog transformation
rlog_counts <- getRlog(counts)

# Save to file
write.table(rlog_counts, "MergedGeneRlog.tsv", sep="\t", quote=FALSE, row.names=FALSE)
