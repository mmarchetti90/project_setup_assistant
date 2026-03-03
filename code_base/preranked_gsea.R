#!/usr/bin/Rscript

# This script reads in differentially expressed genes (DEG) from run_deseq.R and performs enrichment analyses

### ---------------------------------------- ###

parseArgs <- function() {
  
  # Read command line arguments
  args <- commandArgs()
  
  # DEA file location
  dea_path <- args[match("--file_path", args) + 1]
  
  # term2gene path
  term2gene_path <- args[match("--term2gene_path", args) + 1]
  
  # Log2FC threshold
  if("--log2fc_thr" %in% args) {
    
    log2fc <- args[match("--log2fc_thr", args) + 1]
    
  } else {
    
    log2fc <- 0
    
  }
  
  # Adjusted p value threshold
  if("--p_thr" %in% args) {
    
    pval <- args[match("--p_thr", args) + 1]
    
  } else {
    
    pval <- 0.05
    
  }

  # Whether to use all genes or only the significant
  if("--significant_only" %in% args) {

    significant_only <- TRUE

  } else {

    significant_only <- FALSE

  }
  
  return(c(dea_path, term2gene_path, log2fc, pval, significant_only))
  
}

### ---------------------------------------- ###

importDEA <- function(params) {
  
  dea <- read.delim(as.character(params[1]), header = TRUE, row.names = 1, sep = "\t")
  
  dea[is.na(dea$padj), "padj"] <- 1
  
  if(params[5] == TRUE) {
    
    dea <- subset(dea, abs(dea$log2FoldChange) >= as.numeric(params[3]) & dea$padj < as.numeric(params[4]))
    dea["rank"] <- dea$log2FoldChange
    
  } else {
  
    sign <- rep(1, nrow(dea))  
    sign[dea$log2FoldChange < 0] <- -1
    dea["rank"] <- sign * (- log10(dea$padj))
    
  }
  
  return(dea)
  
}

### ---------------------------------------- ###

runGSEA <- function(out, ds, t2g) {
  
  # Creating gene list
  gene_list <- ds$rank
  #names(gene_list) <- as.character(rownames(ds))
  names(gene_list) <- as.character(ds$gene_symbol)
  gene_list <- sort(gene_list, decreasing = TRUE)
  
  # Run analysis
  enrichment <- GSEA(gene_list, pvalueCutoff = 1, pAdjustMethod = "BH", eps = 0, by = "fgsea", minGSSize = 10, maxGSSize = 1000, TERM2GENE = t2g)
  
  if(nrow(enrichment) != 0) {
    
    # Export file
    write.table(enrichment, out, row.names = FALSE, sep='\t')
    
  }
  
  # Plot all terms
  for(term_n in 1:nrow(enrichment)) {

    #if(enrichment@result$p.adjust[term_n] > 0.05) {
    #if(enrichment@result$pvalue[term_n] > 0.05) {
    if(enrichment@result$pvalue[term_n] > 1) {

      next

    }
    
    title_label <- paste(rownames(enrichment@result)[term_n], '\n',
                         'NES = ', enrichment@result$NES[term_n], '\n',
                         'pvalue  = ', enrichment@result$pvalue[term_n], sep='')
                         #'FDR  = ', enrichment@result$p.adjust[term_n], sep='')
    #out_fig <- paste(substring(out, 0, nchar(out) - 4), '_', gsub("/", "-", rownames(enrichment@result)[term_n]), ".png", sep = "")
    #png(file = out_fig, width = 1536, height = 512)
    out_fig <- paste(substring(out, 0, nchar(out) - 4), '_', gsub("/", "-", rownames(enrichment@result)[term_n]), ".pdf", sep = "")
    #pdf(file = out_fig, width = 1536, height = 512)
    gsea_plot <- gseaplot2(enrichment, geneSetID = enrichment$ID[term_n], title = title_label, pvalue_table = FALSE, base_size = 15)
    #print(gsea_plot)
    #dev.off()
    ggsave(file = out_fig, plot = gsea_plot, width = 15, height = 10, units="in", device="pdf", dpi=300)
    
  }
  
  return(enrichment)
  
}

### ------------------MAIN------------------ ###

library(clusterProfiler)
library(enrichplot)
library(RColorBrewer)
library(ggplot2)

parameters <- parseArgs()
print(parameters)

# Import DEA and subset into AllGenes, UpregGenes, and DownregGenes, while filtering for log2fc and padj
dea <- importDEA(parameters)

# Load term2gene matrix for enricher
term2gene <- read.delim(as.character(parameters[2]), header = TRUE, sep = "\t")

# GSEA analysis
if(nrow(dea) != 0) {
  
  gsea_data <- runGSEA("GSEA_AllGenes.tsv", dea, term2gene)
  
}

# Saving datasets to RData object
output_name <- "clusterProfiler_GSEA.RData.gz"
save(parameters, dea, term2gene, gsea_data, file = output_name, compress = T)
