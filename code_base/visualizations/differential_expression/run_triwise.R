#!/usr/bin/Rscript

# This script runs Triwise on gene expression data or on log2 fold change differential expression data

### ---------------------------------------- ###

parseArgs <- function() {
  
  # Read command line arguments
  args <- commandArgs()
  
  # Data matrix
  data_path <- args[match("--data_path", args) + 1]
  
  # Sample names
  sample_names <- unlist(strsplit(args[match("--sample_names", args) + 1], ";"))
  
  # Mode
  if("--expression_mode" %in% args) {
    
    mode <- "expression"
    
  } else if("--fold_change_mode" %in% args) {
    
    mode <- "fold_change"
    
  } else {
    
    mode <- "None"
    
  }
  
  # p value for fold change filtering
  if("--pval_thr" %in% args) {
    
    p_thr <- args[match("--p_thr", args) + 1]
    
  } else {
    
    p_thr <- 0.05
    
  }
  
  # Genes of interest
  if("--goi_path" %in% args) {
    
    goi_path <- args[match("--goi_path", args) + 1]
    goi <- read.delim(goi_path)[, 1]
    
  } else {
    
    goi <- c()
    
  }
  
  # Triwise dotplot rmax
  if("--rmax" %in% args) {
    
    rmax <- args[match("--rmax", args) + 1]
    
  } else {
    
    rmax <- 5
    
  }
  
  return(list(data_path, sample_names, mode, p_thr, goi, rmax))
  
}

### ------------------MAIN------------------ ###

library(triwise)
library(ggplot2)

### Import parameters

parameters <- parseArgs()
print(parameters)
data_path <- parameters[[1]]
sample_names <- parameters[[2]]
mode <- parameters[[3]]
p_thr <- as.numeric(parameters[[4]])
goi <- parameters[[5]]
rmax <- as.numeric(parameters[[6]])

### Load data and get sample names from header

# Load data
data <- read.delim(data_path, sep = "\t", check.names = FALSE)
data[is.na(data)] <- 1

# Remove duplicate gene names/ids and use as row names
data <- data[! duplicated(data[, 1]),]
rownames(data) <- data[, 1]
data <- data[, 2:ncol(data)]

# Format log2 fold changes data
if(mode == "fold_change") {
  
  for(i in seq(1, ncol(data), 2)) {
    
    # Set insignificant fold changes to 1
    data[data[, i+1] >= p_thr, i] <- 1
    
  }
  
  # Remove p value columns
  data <- data[, seq(1, ncol(data), 2)]
  
  # Remove rows that are 0 in all samples
  data <- data[rowSums(data == 0) < 3,]

  # Convert from log2 scale to linear scale
  data <- data**2
  
}

# Rename columns
colnames(data) <- sample_names

### Triwise

# Calculate baricenters
baricenters <- transformBarycentric(data)

# Triwise dotplot and roseplot
if(length(goi) == 0) {
  
  # Dotplot
  plotDotplot(baricenters,
              Coi = sample_names,
              rmax = rmax) +
    ylim(- rmax - .25, rmax + .25) +
    xlim(- rmax - .25, rmax + .25)
  ggsave(filename="triwise_dotplot.png", bg = "white", dpi=300)
  
  # Roseplot
  plotRoseplot(baricenters,
               Coi = sample_names,
               rmax = 5)
  ggsave(filename="triwise_roseplot.png", bg = "white", dpi=300)
  
} else {
  
  # Dotplot
  plotDotplot(baricenters,
              Coi = sample_names,
              Goi = goi,
              rmax = rmax) +
    ylim(- rmax - .25, rmax + .25) +
    xlim(- rmax - .25, rmax + .25)
  ggsave(filename="triwise_dotplot.png", bg = "white", dpi=300)
  
  # Roseplot
  plotRoseplot(baricenters,
               Coi = sample_names,
               Goi = goi)
  ggsave(filename="triwise_roseplot.png", bg = "white", dpi=300)
  
}
