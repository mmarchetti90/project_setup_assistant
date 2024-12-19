#!/usr/bin/Rscript

# This script runs MuSiC deconvolution
# Can use MuSiC 1 or 2 (MuSiC 2 never worked in my hands)

### ---------------------------------------- ###

parseArgs <- function() {
  
  # Read command line arguments
  args <- commandArgs()
  
  # Paths to data
  bulk_counts_file <- args[match("--bulk_counts", args) + 1]
  bulk_metadata_file <- args[match("--bulk_metadata", args) + 1]
  single_cell_counts_file <- args[match("--single_cell_counts", args) + 1]
  single_cell_metadata_file <- args[match("--single_cell_metadata", args) + 1]
  
  # Bulk control condition for MuSiC2 TOAST
  if("--bulk_ctrl_condition" %in% args) {
    
    bulk_ctrl_condition <- args[match("--bulk_ctrl_condition", args) + 1]
    
  } else {
    
    bulk_ctrl_condition <- NA
    
  }
  
  # Run music_prop (music_version=1) or music2_prop (music_version=2) or music2_prop_toast (music_version=toast)
  music_version <- args[match("--music_version", args) + 1]
  
  return(c(bulk_counts_file, bulk_metadata_file, single_cell_counts_file, single_cell_metadata_file, bulk_ctrl_condition, music_version))
  
}

### ---------------------------------------- ###

createExpressionSet <- function(expression_file, metadata_file) {
  
  expression <- as.matrix(read.delim(expression_file, header=TRUE, row.names=1, sep='\t'))
  
  metadata <- read.delim(metadata_file, header=TRUE, row.names=1, sep='\t')
  metadata <- new("AnnotatedDataFrame", data=metadata)
  
  expression_set <- ExpressionSet(assayData=expression, phenoData=metadata)
  
  return(expression_set)
  
}

### ---------------------------------------- ###

createSingleCellExperiment <- function(expression_file, metadata_file) {
  
  expression <- as.matrix(read.delim(expression_file, header=TRUE, row.names=1, sep='\t', check.names=FALSE))
  
  metadata <- read.delim(metadata_file, header=TRUE, row.names=1, sep='\t')
  
  sce <- SingleCellExperiment(list(counts=expression),
                              colData=metadata)
  
  
  return(sce)
  
}

### ---------------------------------------- ###

runMuSiC <- function(bulk_expr, bulk_metadata, bulk_ctrl_condition, single_cell_data, music_version) {
  
  set.seed(42)
  
  if(music_version == "1") {
    
    print("Running MuSiC 1")
    
    estimation <- music_prop(
      bulk_expr,
      single_cell_data,
      markers = NULL,
      clusters = 'cellType',
      samples = 'sampleID',
      select.ct = unique(single_cell_data@colData@listData$cellType),
      cell_size = NULL,
      ct.cov = FALSE,
      verbose = TRUE,
      iter.max = 1000,
      nu = 1e-04,
      eps = 0.01,
      centered = FALSE,
      normalize = FALSE
    )
    
  } else if(music_version == "2") {
    
    print("Running MuSiC 2")
    
    estimation <- music2_prop(
      bulk_expr[, colnames(bulk_expr) %in% rownames(bulk_metadata[bulk_metadata$genotype == parameters[5],])],
      bulk_expr[, colnames(bulk_expr) %in% rownames(bulk_metadata[bulk_metadata$genotype != parameters[5],])],
      sc.sce = single_cell_data,
      clusters = 'cellType',
      samples = 'sampleID',
      select.ct = unique(single_cell_data@colData@listData$cellType),
      expr_low = 20,
      prop_r = 0.1,
      eps_c = 0.05,
      eps_r = 0.01,
      n_resample = 20,
      sample_prop = 0.5,
      cutoff_expr = 0.05,
      cutoff_c = 0.05,
      cutoff_r = 0.01,
      maxiter = 200,
      markers = NULL,
      cell_size = NULL,
      ct.cov = FALSE,
      centered = FALSE,
      normalize = FALSE
    )
    
  } else if(music_version == "toast") {
    
    print("Running MuSiC 2 TOAST")
    
    estimation <- music2_prop_toast(
      bulk_expr[, colnames(bulk_expr) %in% rownames(bulk_metadata[bulk_metadata$genotype == parameters[5],])],
      bulk_expr[, colnames(bulk_expr) %in% rownames(bulk_metadata[bulk_metadata$genotype != parameters[5],])],
      sc.sce = single_cell_data,
      clusters = 'cellType',
      samples = 'sampleID',
      select.ct = unique(single_cell_data@colData@listData$cellType),
      expr_low = 20,
      prop_r = 0.1,
      eps_c = 0.05,
      eps_r = 0.01,
      cutoff_c = 10^(-3),
      cutoff_r = 10^(-3),
      cap = 0.3,
      maxiter = 200
    )
    
  } else {
    
    print("ERROR: invalid MuSiC version specified")
    estimation <- NA
    
  }
  
  return(estimation)
  
}

### ---------------------------------------- ###

plot_estimates_heatmap <- function(prop, row_annot, outname) {
  
  # Sort cell types by name
  prop <- prop[, order(colnames(prop))]
  
  # Create color palette
  heatmap_colors <- colorRampPalette(c("white", "red"))(50)
  
  # Heatmap
  pheatmap(prop,
           border_color='black',
           color=heatmap_colors,
           annotation_row=row_annot,
           cellwidth=15,
           cellheight=15,
           scale="none",
           cluster_rows=TRUE,
           cluster_cols=FALSE,
           legend=TRUE,
           angle_col=45,
           fontsize=12,
           filename=outname)
  
}

### ---------------------------------------- ###

plot_estimates_barplot <- function(prop, metadata, outname_prefix) {
  
  for(i in 1:ncol(prop)) {
    
    plot_data <- cbind(metadata, prop[rownames(metadata), i])
    colnames(plot_data) <- c('sample_id', 'genotype', 'proportion')
    
    ggplot(data=plot_data, aes(x=sample_id, y=proportion, fill=genotype)) +
      geom_bar(stat="identity", color="black") +
      xlab(NULL) +
      ylab("Population proportion (%)") +
      scale_x_discrete(limits=factor(plot_data$sample_id), labels=plot_data$genotype) +
      theme(axis.text.x=element_text(angle=45, vjust=1, hjust=1))
    ggsave(outname <- gsub("/", "-", paste(outname_prefix, colnames(prop)[i], ".pdf", sep="")))
    
  }
  
}

### ---------------------------------------- ###

check_genotype_effect <- function(prop, metadata, bulk_ctrl_condition) {
  
  # Create a linear model to evaluate the effect of genotype on populations frequency
  
  # Define p value threshold
  p_thr <- 0.05
  
  # Init flags
  flags <- c()
  
  for(cell_type in colnames(prop)) {
    
    # Prepare data
    data <- as.data.frame(cbind(metadata$genotype, prop[rownames(metadata), cell_type]))
    colnames(data) <- c("genotype", "proportion")
    data$genotype <- relevel(as.factor(data$genotype), ref=bulk_ctrl_condition) # This makes it so genotype effects are computed in respect ot control
    
    # Compute model
    model <- lm(proportion ~ genotype, data)
    
    # Summarize model
    coefficients <- summary(model)$coefficients
    rownames(coefficients) <- gsub("genotype", "", rownames(coefficients))
    colnames(coefficients) <- c("LinearCoefficient", "StdErr", "t-value", "p-value")
    
    # Write to file
    outname <- paste(cell_type, "_genotype_effect.tsv", sep="")
    outname <- gsub(" ", "_", outname)
    outname <- gsub("/", "-", outname)
    write.table(coefficients, outname, sep="\t", quote=FALSE)
    
    # Check for significant enrichment/depletions
    for(significant_effects in which(coefficients[, 4] < p_thr)) {
      
      genotype <- rownames(coefficients)[significant_effects]
      
      if(genotype == "(Intercept)" | is.null(significant_effects)) {
        
        next
        
      }
      
      base_frequency <- median(data[data$genotype == bulk_ctrl_condition, "proportion"])
      genotype_frequency <- median(data[data$genotype == genotype, "proportion"])
      
      if(genotype_frequency > base_frequency) {
        
        new_flag <- paste(cell_type, " enriched in ", genotype, " (", genotype_frequency, "%) vs ", bulk_ctrl_condition, " (", base_frequency, "%)", sep="")
        
      } else {
        
        new_flag <- paste(cell_type, " depleted in ", genotype, " (", genotype_frequency, "%) vs ", bulk_ctrl_condition, " (", base_frequency, "%)", sep="")
        
      }
      
      flags <- c(flags, new_flag)
      
    }
    
  }
  
  # Save flags to file
  flags_file <- file("cell_type_flags.txt")
  writeLines(flags, flags_file)
  close(flags_file)
  
}

### ------------------MAIN------------------ ###

library(Biobase)
library(MuSiC)
library(pheatmap)
library(SingleCellExperiment)

### LOAD DATA ------------------------------ ###

# Reading parameters
parameters <- parseArgs()

# Load bulk counts and metadata
bulk_expr <- as.matrix(read.delim(parameters[1], header=TRUE, row.names=1, sep='\t'))
bulk_metadata <- read.delim(parameters[2], header=TRUE, row.names=1, sep='\t')

# Load single-cell counts and metadata and convert to ExpressionSet
single_cell_data <- createSingleCellExperiment(parameters[3], parameters[4])

### ESTIMATE CELL PROPORTIONS -------------- ###

# Run MuSiC
music_estimation <- runMuSiC(bulk_expr, bulk_metadata, parameters[5], single_cell_data, parameters[6])

# Save as RDS object
saveRDS(music_estimation, "MuSiC_Estimation.rds")

# Export proportion matrices
write.table(music_estimation$Est.prop.weighted, "estimated_proportions.tsv", quote=FALSE, sep='\t')

### PLOTTING ------------------------------- ###

# Prepare heatmap sample annotation
sample_annotation <- data.frame(genotype = bulk_metadata$genotype, row.names = rownames(bulk_metadata))
sample_annotation$genotype <- factor(sample_annotation$genotype, levels = unique(sample_annotation$genotype))

# Plotting heatmap
plot_estimates_heatmap(music_estimation$Est.prop.weighted, sample_annotation, 'estimated_proportions_heatmap.pdf')

# Bar plots
plot_estimates_barplot(music_estimation$Est.prop.weighted, bulk_metadata, "estimated_proportions_barplot_")

### GENOTYPE EFFECT ------------------------ ###

if(! is.na(parameters[5])) {
  
  check_genotype_effect(music_estimation$Est.prop.weighted, bulk_metadata, parameters[5])
  
}
