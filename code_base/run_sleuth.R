#!/usr/bin/Rscript

# This script reads in Kallisto transcript counts and performs a differential splicing analysis
# N.B. This is designed for simple binary comparisons as described in a comparison design file in the column "condition" and field "reference"

### ---------------------------------------- ###

parseArgs <- function() {
  
  # Read command line arguments
  args <- commandArgs()
  
  # Path to Kallisto output directories
  kallisto_out_path <- args[match("--kallisto_out_path", args) + 1]
  
  # Design file
  design_file <- args[match("--design", args) + 1]
  
  # Cores to use for sleuth_prep
  if("--num_cores" %in% args) {
    
    num_cores <- args[match("--num_cores", args) + 1]
    
  } else {
    
    num_cores <- 1
    
  }

  return(c(kallisto_out_path, design_file, num_cores))
  
}

### ---------------------------------------- ###

importDesign <- function(design_file, kallisto_dirs, kallisto_dirs_path) {
  
  # Reading design file
  temp <- read.delim(design_file, header = FALSE, sep = "\t", stringsAsFactors = FALSE)
  
  # Extracting analysis name
  analysis_name <- as.character(temp[1, 2])
  
  # Extracting comparison formula
  condition <- as.formula(temp[2, 2])
  
  # Extracting reference
  reference <- as.character(temp[3, 2])
  
  # Generating sample_info dataframe
  sample_info <- as.data.frame(temp[5 : nrow(temp), 1 : ncol(temp)], stringsAsFactors = TRUE)
  colnames(sample_info) <- temp[4, 1 : ncol(temp)]
  sample_info <- addKallistoPaths(kallisto_dirs, kallisto_dirs_path, sample_info)
  
  # Releveling "condition" to set reference
  sample_info$condition <- relevel(as.factor(sample_info$condition), ref=reference)
  
  return(list("analysis_name" = analysis_name, "formula" = condition, "reference" = reference, "sample_info" = sample_info))
  
}

### ---------------------------------------- ###

addKallistoPaths <- function(dirs, dirs_path, data) {
  
  paths <- c()
  for(sample in data$sample) {
    
    if(paste(dirs_path, "/", sample, "_kallisto_out", sep="") %in% dirs) {
      
      new_path <- paste(dirs_path, "/", sample, "_kallisto_out", sep="")
      
    } else {
      
      new_path <- "NA"
      
    }
    
    paths <- c(paths, new_path)
    
  }
  
  data["path"] <- paths
  data <- data[data$path != "NA",]
  
  data <- data[order(data$condition, data$sample),]
  
  return(data)
  
}

### ---------------------------------------- ###

runDEA <- function(params, des, cores) {
  
  # Create Sleuth object
  sleuth_obj <- sleuth_prep(des$sample_info, extra_bootstrap_summary = TRUE, num_cores = cores)
  
  # Fit Sleuth models (full and reduced)
  sleuth_obj <- sleuth_fit(sleuth_obj, des$formula, 'full') # Full model
  sleuth_obj <- sleuth_fit(sleuth_obj, ~1, 'reduced') # Intercept only model
  
  # Differential transcript expression testing using Likelihood Ratio Test
  sleuth_obj <- sleuth_lrt(sleuth_obj, 'reduced', 'full')
  
  # Differential transcript expression testing using Wald Test (to get a fold change, which lrt does not provide)
  beta_name <- paste("condition", des$sample_info[des$sample_info$condition != des$reference, "condition"][1], sep="")
  sleuth_obj <- sleuth_wt(sleuth_obj, beta_name, "full")
  
  # Extracting results
  results_lrt <- sleuth_results(sleuth_obj, "reduced:full", test_type = "lrt", show_all = FALSE)
  results_wd <- sleuth_results(sleuth_obj, beta_name, test_type = 'wt', show_all = FALSE)
  
  # Saving results and sleuth object
  write.table(results_lrt, paste("Sleuth_LRT_", des$analysis_name, ".tsv", sep=""), row.names = FALSE, sep = "\t")
  write.table(results_wd, paste("Sleuth_Wald_", des$analysis_name, ".tsv", sep=""), row.names = FALSE, sep = "\t")
  saveRDS(sleuth_obj, paste("Sleuth_", des$analysis_name, ".rds", sep = ""))
  
  ### QUALITY CONTROL PLOTS ---------------- ###
  
  # Group density plot
  output_name <- paste("Sleuth_", des$analysis_name, "_GroupDensity.png", sep = "")
  plot_group_density(sleuth_obj, use_filtered = TRUE,
                     units = "est_counts",
                     trans = "log",
                     grouping = setdiff(colnames(sleuth_obj$sample_to_covariates), "sample"),
                     offset = 1)
  ggsave(filename=output_name, width = 20, height = 15, units = "cm", dpi=300)
  
  # MA plot WT (LRT not supported)
  output_name <- paste("Sleuth_Wald", des$analysis_name, "_MA-Plot.png", sep = "")
  plot_ma(sleuth_obj, 
          beta_name, 
          test_type = "wt", 
          which_model = "full",
          sig_level = 0.05, 
          point_alpha = 0.2, 
          sig_color = "red",
          highlight = NULL, 
          highlight_color = "green")
  ggsave(filename=output_name, width = 20, height = 15, units = "cm", dpi=300)
  
  # Mean-variance relationship plot
  output_name <- paste("Sleuth_", des$analysis_name, "_MeanVsVariance.png", sep = "")
  plot_mean_var(sleuth_obj,
                which_model = "full",
                point_alpha = 0.4,
                point_size = 2,
                point_colors = c("black", "dodgerblue"),
                smooth_alpha = 1,
                smooth_size = 0.75,
                smooth_color = "red")
  ggsave(filename=output_name, width = 20, height = 15, units = "cm", dpi=300)
  
  # PCA plot
  output_name <- paste("Sleuth_", des$analysis_name, "_PCA.png", sep = "")
  plot_pca(sleuth_obj,
           pc_x = 1L,
           pc_y = 2L,
           use_filtered = TRUE,
           units = "est_counts",
           text_labels = TRUE,
           color_by = des$sample_info$condition,
           point_size = 3,
           point_alpha = 0.8)
  ggsave(filename=output_name, width = 30, height = 30, units = "cm", dpi=300)
  
  # Scree plot
  output_name <- paste("Sleuth_", des$analysis_name, "_ScreePlot.png", sep = "")
  plot_pc_variance(sleuth_obj,
                   use_filtered = TRUE,
                   units = "est_counts",
                   pca_number = NULL,
                   scale = FALSE,
                   PC_relative = NULL)
  ggsave(filename=output_name, width = 20, height = 15, units = "cm", dpi=300)
  
  # QQ plot for LRT
  output_name <- paste("Sleuth_LRT_", des$analysis_name, "_QQ.png", sep = "")
  plot_qq(sleuth_obj,
          "reduced:full",
          test_type = "lrt",
          which_model = "full",
          sig_level = 0.05,
          point_alpha = 0.2,
          sig_color = "red",
          highlight = NULL,
          highlight_color = "green",
          line_color = "blue")
  ggsave(filename=output_name, width = 20, height = 15, units = "cm", dpi=300)
  
  # QQ plot for WT
  output_name <- paste("Sleuth_Wald_", des$analysis_name, "_QQ.png", sep = "")
  plot_qq(sleuth_obj,
          beta_name,
          test_type = "wt",
          which_model = "full",
          sig_level = 0.05,
          point_alpha = 0.2,
          sig_color = "red",
          highlight = NULL,
          highlight_color = "green",
          line_color = "blue")
  ggsave(filename=output_name, width = 20, height = 15, units = "cm", dpi=300)
  
  # Heatmap of top 50 transcripts
  transcripts_to_plot <- results_lrt$target_id[1:min(50, nrow(results_lrt))]
  
  if(length(transcripts_to_plot) > 10) {
   
    output_name <- paste("Sleuth_", design$analysis_name, "_Heatmap.png", sep = "")
    png(file = output_name, width = 800, height = 600)
    plot_transcript_heatmap(sleuth_obj,
                            transcripts=transcripts_to_plot,
                            units = "tpm",
                            trans = "log",
                            cluster_transcripts = TRUE,
                            offset = 1,
                            color_high = "#581845",
                            color_mid = "#FFC300",
                            color_low = "#DAF7A6",
                            x_axis_angle = 50,
                            annotation_cols = setdiff(colnames(sleuth_obj$sample_to_covariates), "sample"))
    dev.off()
    
  }
  
  # Volcano plot for Wald Test (LRT not supported)
  output_name <- paste("Sleuth_Wald_", des$analysis_name, "_VolcanoPlot.png", sep = "")
  plot_volcano(sleuth_obj,
               beta_name,
               test_type = "wt",
               which_model = "full",
               sig_level = 0.05,
               point_alpha = 0.2,
               sig_color = "red",
               highlight = NULL)
  ggsave(filename=output_name, width = 20, height = 15, units = "cm", dpi=300)
  
}

### ------------------MAIN------------------ ###

library(sleuth)
library(ggplot2)

parameters <- parseArgs()

# Import counts and experimental design
kallisto_outputs <- list.files(path=as.character(parameters[1]), pattern="*_kallisto_out", full.names = TRUE, include.dirs = TRUE)
design <- importDesign(as.character(parameters[2]), kallisto_outputs, parameters[1])

# Running DEA
runDEA(parameters, design, as.integer(parameters[3]))
