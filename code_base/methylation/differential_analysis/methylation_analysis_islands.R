#!/usr/bin/Rscript

# This script analyzes methylation coverage files from Bismark at a base-pair or CpG-island level
# N.B. Current version only works for hg38

### ---------------------------------------- ###

parseArgs <- function() {
  
  # Read command line arguments
  args <- commandArgs()
  
  # Path to counts file directory
  path_to_files <- args[match("--cov_dir", args) + 1]
  
  # Path to design file
  if("--comparison_file" %in% args) {
    
    samples_file <- args[match("--comparison_file", args) + 1]
    
  } else {
    
    samples_file <- ""
    
  }
  
  # CPUs
  if("--cpus" %in% args) {
    
    cpus <- as.numeric(args[match("--cpus", args) + 1])
    
  } else {
    
    cpus <- 1
    
  }
  
  # Significance threshold
  if("--significance_threshold" %in% args) {
    
    significance_threshold <- as.numeric(args[match("--significance_threshold", args) + 1])
    
  } else {
    
    significance_threshold <- 0.05
    
  }
  
  # Significance threshold
  if("--use_pvalue" %in% args) {
    
    significance_parameter <- 'pvalue'
    
  } else {
    
    significance_parameter <- 'padj'
    
  }

  # Significance threshold
  if("--island_analysis" %in% args) {
    
    analysis_level <- 'island'
    
  } else {
    
    analysis_level <- 'bp'
    
  }
  
  # Methylation analysis type
  if("--utest" %in% args) {
    
    analysis_type <- 'utest'
    
  } else {
    
    analysis_type <- 'methylkit'
    
  }
  
  # Paired samples toggle
  if("--paired" %in% args) {
    
    paired <- TRUE
    
  } else {
    
    paired <- FALSE
    
  }
  
  # Islands bp extension upstream and downstream prior to annotation
  if("--island_extension" %in% args) {
    
    island_extension <- as.numeric(args[match("--island_extension", args) + 1])
    
  } else {
    
    island_extension <- 0
    
  }
  
  # Islands bp extension upstream and downstream prior to annotation
  if("--correct_batch" %in% args) {
    
    correct_batch <- TRUE
    
  } else {
    
    correct_batch <- FALSE
    
  }
  
  return(c(path_to_files, samples_file, cpus, significance_threshold, significance_parameter, analysis_level, analysis_type, paired, island_extension, correct_batch))
  
}

### ---------------------------------------- ###

structureWorkDir <- function(samples_file, analysis_level) {
  
  # Defining basedir name
  if(samples_file == "") {

    if(analysis_level == "island") {

      basedir <- "CpG_Islands_Analysis"

    } else {

      basedir <- "DMC_Analysis"

    }
    
  } else {

    if(analysis_level == "island") {

      suffix <- "_CpG_Islands_Analysis"

    } else {

      suffix <- "_DMC_Analysis"

    }
    
    basedir <- substr(basename(samples_file),
                      regexpr("Comparison_", basename(samples_file)) + 11,
                      regexpr(".tsv", basename(samples_file)) - 1)
    basedir <- paste(basedir, suffix, sep = "")
    
  }
  
  for(subdir in c("RelationshipsPlot", "DifferentialMethylation_Plots", "DifferentialMethylation_Analysis", "MethylDB")) {
      
      dir.create(paste(basedir, subdir, sep='/'), recursive = T, showWarnings = F)
    
  }

  return(basedir)
  
}

### ---------------------------------------- ###

getSamplesInfo <- function(path_to_files, samples_file, cov_files) {
  
  if(samples_file == "") { # samples file was not provided
    
    info <- DataFrame(SampleID = as.vector(sapply(cov_files, FUN = function(x) { substr(basename(x), 0, regexpr("_", basename(x))[1] - 1) } )),
                      FilePath = as.vector(sapply(cov_files, FUN = function(x) { paste(path_to_files, x, sep = "/") } )),
                      Condition = factor(rep(0, length(cov_files)), levels=c(0, 1)))
    
  } else {
    
    info <- read.delim(samples_file, sep = "\t", stringsAsFactors = F)
    info$Condition <- as.integer(as.factor(info$Condition)) - 1
    
  }
  
  return(info)
  
}

### ---------------------------------------- ###

parseCoverageFiles <- function(basedir, info) {
  
  methylation_data <- methRead(as.list(info$FilePath),
                               sample.id = as.list(info$SampleID),
                               assembly = "hg18",
                               treatment = info$Condition,
                               header = FALSE,
                               context = "CpG",
                               mincov = 10,
                               dbtype = "tabix",
                               dbdir = paste(basedir, "MethylDB", sep='/'),
                               pipeline = "bismarkCoverage")
  
  return(methylation_data)
  
}

### ---------------------------------------- ###

plotSamplesRelationship <- function(basedir, methylation_data) {
  
  # Plot sample correlation
  #png(filename = paste(basedir, "/RelationshipsPlot/SamplesCorrelation.png", sep=""), width = 1024, height = 1024)
  #getCorrelation(methylation_data, plot = T)
  #dev.off()
  
  # Clustering
  png(filename = paste(basedir, "/RelationshipsPlot/SamplesClusters.png", sep=""), width = 512, height = 512)
  clusterSamples(methylation_data, dist = "correlation", method = "ward.D2", plot = T)
  dev.off()
  
  # PCA Screeplot
  png(filename = paste(basedir, "/RelationshipsPlot/SamplesPCAScreeplot.png", sep=""), width = 512, height = 512)
  PCASamples(methylation_data, screeplot = T)
  dev.off()
  
  # PCA 1 vs 2
  png(filename = paste(basedir, "/RelationshipsPlot/SamplesPCA.png", sep=""), width = 512, height = 512)
  pca_data = PCASamples(methylation_data, obj.return=TRUE)
  dev.off()
  
  # Saving PCA data to file
  write.table(pca_data$x, paste(basedir, "/RelationshipsPlot/pca_data.tsv", sep=""), col.names=NA, sep="\t")
  
}

### ---------------------------------------- ###

batch_correction <- function(si, md, bd, var, params_toggle) {
  
  # Check if batch correction is needed
  if(var %in% colnames(si)) {
    
    if(length(unique(si[[var]])) > 1) {
      
      for(b in unique(si[[var]])) {
        
        if(sum(si[[var]] == b) > 1) {
          
          batch_correction_toggle <- TRUE
          
        } else {
          
          batch_correction_toggle <- FALSE
          break
          
        }
        
      }
      
    } else {
      
      batch_correction_toggle <- FALSE
      
    }
    
  } else {
    
    batch_correction_toggle <- FALSE
    
  }
  
  if(params_toggle == FALSE) {
    
    batch_correction_toggle <- FALSE
    
  }
  
  # Correcting for batch effects, using methylKit
  # Association between PC and a variables are used to filter out those PCs that associate above a significance threshold,
  # then data is reconstructed using only the PCs that passed filtering
  if((batch_correction_toggle) & (var %in% colnames(si))) {
    
    print(paste("### Correcting", var, "effect", sep=" "))
    
    # Gathering batch info
    batch_annotation <- data.frame('Batch' = si[match(md@sample.ids, si$SampleID), var])
    
    # Computing association of Batch with principal components
    association <- assocComp(md, batch_annotation)
    
    # Saving results to file
    write.table(association$vars, paste(bd, '/RelationshipsPlot/pca_variance_', var, '.tsv', sep=''), row.names=F, col.names=F, sep='\t', quote=FALSE)
    write.table(association$association, paste(bd, '/RelationshipsPlot/pca_', var, '_association_pval.tsv', sep=''), row.names=F, col.names=F, sep='\t', quote=FALSE)
    
    # Removing components associated with batch
    md <- removeComp(md, comp=which(association$association < 0.05))
    
  }
  
  return(md)
  
}

### ---------------------------------------- ###

plotDMC <- function(basedir, dmc_analysis) {
  
  # Basic volcano plot
  png(filename = paste(basedir, "/DifferentialMethylation_Plots/VolcanoPlot.png", sep=""), width = 1024, height = 1024)
  plot(dmc_analysis$meth.diff, -log10(dmc_analysis$qvalue), type="p", pch=21, col="black", bg="gray",
       main=paste(basedir, "_Volcano", sep=""), xlab="Methylation difference", ylab="- Log10(qvalue)")
  dev.off()
  
  # Frequency of hyper- and hypo-methylation events per chromosome
  png(filename = paste(basedir, "/DifferentialMethylation_Plots/DifferentialMethylationFrequency.png", sep=""), width = 1024, height = 512)
  diffMethPerChr(dmc_analysis)
  dev.off()
  
}

### ---------------------------------------- ###

extractDMC <- function(basedir, dmc_analysis, significance_threshold, significance_parameter) {
  
  # Extracting differentially methylated regions
  if(significance_parameter == 'pvalue') {
    
    dmc_regions <- dmc_analysis[dmc_analysis$pvalue <= significance_threshold]
    
  } else {
    
    dmc_regions <- dmc_analysis[dmc_analysis$qvalue <= significance_threshold]
    
  }
  
  # Saving to file
  write.table(dmc_regions, paste(basedir, "/DifferentialMethylation_Analysis/DMC_Regions.tsv", sep=""), row.names=F, sep="\t")
  
  return(dmc_regions)
  
}

### ---------------------------------------- ###

plotMethylationTSS <- function(basedir, dmc_regions, extension=5000, palette_clipping=50) {
  
  # Extracting genes and their tss as Granges objects
  genes <- genes(TxDb.Hsapiens.UCSC.hg38.knownGene)
  tss <- promoters(genes, upstream = 0, downstream = 1)
  
  # Adding "chr" in front of every dmc_regions chromosome name to match the grange in tss
  dmc_regions$chr <- paste("chr", dmc_regions$chr, sep="")
  
  # Converting dmc_regions to GRanges for annotatr
  dmc_regions_gr <- as(dmc_regions, "GRanges")
  
  # Adding an abs.meth.diff column to dmc_regions_gr
  dmc_regions_gr$abs.meth.diff <- abs(dmc_regions_gr$meth.diff)
  
  # Defining a color palette for meth.diff and abs.meth.diff
  if(min(dmc_regions_gr$meth.diff) > 0) {
    
    palette1 <- colorRamp2(c(0, min(max(dmc_regions_gr$meth.diff), palette_clipping)), c("white", "red"))
    
  } else {
    
    palette1 <- colorRamp2(c(max(min(dmc_regions_gr$meth.diff), -palette_clipping), 0, min(max(dmc_regions_gr$meth.diff), palette_clipping)), c("blue", "white", "red"))
    
  }
  palette2 <- colorRamp2(c(0, min(max(dmc_regions_gr$abs.meth.diff), palette_clipping)), c("white", "red"))
  
  # Plotting heatmap of methylation differences distribution around TSS
  normalized_matrix <- normalizeToMatrix(dmc_regions_gr, tss, value_column = "meth.diff", extend = extension, background = 0, mean_mode="absolute", w=50)
  
  png(filename = paste(basedir, "/DifferentialMethylation_Plots/TSSMethylationDifference_AllGenes.png", sep=""), width = 768, height = 1536, res=150)
  plot <- EnrichedHeatmap(normalized_matrix, col = palette1, name = "Methylation difference")
  print(plot)
  dev.off()
  
  # Removing genes with no DMCs and replotting
  normalized_matrix <- normalized_matrix[rowSums(as.data.frame(normalized_matrix)) != 0,]
  
  png(filename = paste(basedir, "/DifferentialMethylation_Plots/TSSMethylationDifference_WithDMCOnly.png", sep=""), width = 768, height = 1536, res=150)
  plot <- EnrichedHeatmap(normalized_matrix, col = palette1, name = "Methylation difference")
  print(plot)
  dev.off()
  
  # Plotting heatmap of absolute methylation differences distribution around TSS
  normalized_matrix <- normalizeToMatrix(dmc_regions_gr, tss, value_column = "abs.meth.diff", extend = extension, background = 0, mean_mode="absolute", w=50)
  
  png(filename = paste(basedir, "/DifferentialMethylation_Plots/TSSMethylationDifference_Absolute_All_genes.png", sep=""), width = 768, height = 1536, res=150)
  plot <- EnrichedHeatmap(normalized_matrix, col = palette2, name = "Methylation difference")
  print(plot)
  dev.off()
  
  # Removing genes with no DMCs and replotting
  normalized_matrix <- normalized_matrix[rowSums(as.data.frame(normalized_matrix)) != 0,]
  
  png(filename = paste(basedir, "/DifferentialMethylation_Plots/TSSMethylationDifference_Absolute_WithDMCOnly.png", sep=""), width = 768, height = 1536, res=150)
  plot <- EnrichedHeatmap(normalized_matrix, col = palette2, name = "Methylation difference")
  print(plot)
  dev.off()
  
}

### ---------------------------------------- ###

annotateDMC <- function(basedir, dmc_regions, island_extension) {
  
  # Adding "chr" in front of every dmc_regions chromosome name to match the grange created by build_annotations
  dmc_regions$chr <- paste("chr", dmc_regions$chr, sep="")
  
  # Extending islands
  dmc_regions$start <- dmc_regions$start - island_extension
  dmc_regions[dmc_regions$start < 0, 'start'] <- 0
  dmc_regions$end <- dmc_regions$end + island_extension
  
  # Converting dmc_regions to GRanges for annotatr
  dmc_regions_gr <- as(dmc_regions, "GRanges")
  
  # Building annotations
  annotations = build_annotations(genome = 'hg38', annotations = c('hg38_cpgs', 'hg38_basicgenes'))
  
  # Annotating DMC regions
  dmc_annotation <- annotate_regions(regions = dmc_regions_gr, annotations = annotations, ignore.strand = T, quiet = F)
  
  # Saving to file
  write.table(data.frame(dmc_annotation), paste(basedir, "/DifferentialMethylation_Analysis/DMC_Annotation.tsv", sep=""), row.names=F, sep="\t")
  
  # Plotting a summary of the annotations
  #annotation_summary <- summarize_annotations(annotated_regions = dmc_annotation, quiet = TRUE)
  summary_plot <- plot_annotation(annotated_regions = dmc_annotation,
                                  plot_title = 'Annotation of differentially methylated cytosines',
                                  x_label = '',
                                  y_label = 'Count')
  ggsave(paste(basedir, "/DifferentialMethylation_Plots/AnnotationSummary.png", sep=""), plot = summary_plot, dpi=300)
  
}

### ------------------MAIN------------------ ###

library(methylKit)
library(graphics)
library(annotatr)
library(AnnotationHub)
library(TxDb.Hsapiens.UCSC.hg38.knownGene)
library(GenomicFeatures)
library(EnrichedHeatmap)
library(ggplot2)
library(circlize)
library(stringr)
library(sva)
library(matrixStats)

### LOAD INFO AND STRUCTURE RESULTS DIR ---- ###

print("### Reading info")

# Reading parameters and listing coverage files
parameters <- parseArgs()
print(parameters)
cov_files <- list.files(path=parameters[1], pattern="*.cov*", full.names = F)

# Structuring workdir
basedir <- structureWorkDir(parameters[2], parameters[6])

### LOAD DATA ----------------------------- ###

# Import samples info
print("### Import samples' info")
samples_info <- getSamplesInfo(parameters[1], parameters[2], cov_files)

# Import methylation data
print("### Reading sample data into methylRawListDB object")
methylation_data <- parseCoverageFiles(basedir, samples_info)

# Save raw data
saveRDS(methylation_data, paste(basedir, '/MethylDB/raw_data.rds', sep=""))

# Merging cytosine counts for island-level analysis (if option is selected)
if(parameters[6] == "island") {

  # Create GRanges object of CpG islands
  print("### Creating CpG island GRange object")
  cpg_islands_coordinates <- build_annotations(genome = 'hg38', annotations = 'hg38_cpgs')

  # Only select islands (no shores, or other)
  cpg_islands_coordinates <- cpg_islands_coordinates[cpg_islands_coordinates$type == "hg38_cpg_islands",]
  
  # Converting Granges object chromosome names to match the methylRawList
  cpg_islands_coordinates <- as.data.frame(cpg_islands_coordinates)
  cpg_islands_coordinates$seqnames <- str_replace(cpg_islands_coordinates$seqnames, 'chr', '')
  cpg_islands_coordinates <- as(cpg_islands_coordinates, "GRanges")
  
  # N.B. For whatever reason, the function fails if run with save.db = TRUE
  print("### Merging cytosine counts to CpG-islands")
  #methylation_data <- regionCounts(methylation_data, cpg_islands_coordinates, save.db = TRUE, suffix = "island_level")
  methylation_data <- regionCounts(methylation_data, cpg_islands_coordinates, save.db = FALSE)
  
  # Save raw data after merging
  saveRDS(methylation_data, paste(basedir, '/MethylDB/raw_data_merged.rds', sep=""))

}

# Filtering too low/high coverage (risk of PCR bias and poor statistics, respectively)
print("### Filtering cytosines with too low/high coverage")
methylation_data <- filterByCoverage(methylation_data, lo.count = 10, lo.perc = NULL, hi.count = NULL, hi.perc = 99.9)

# Normalization
print("### Data normalization")
methylation_data <- normalizeCoverage(methylation_data, method = "median")

# Unite sample data for downstream differential analyses (no need to destrand)
print("### Uniting data")
methylation_data <- unite(methylation_data, destrand = FALSE)

# Save object prior to batch correction
saveRDS(methylation_data, paste(basedir, '/MethylDB/united_data.rds', sep=""))

# Correction for batch and gender effects
for(var in c("Batch", "Gender")) {
  
  methylation_data <- batch_correction(samples_info, methylation_data, basedir, var, parameters[10])
  
}

# Save object after batch correction
saveRDS(methylation_data, paste(basedir, '/MethylDB/united_batch_corrected_data.rds', sep=""))

# Plotting samples correlation, clustering, and PCA 1,2
print("### Plotting sample relationships")
plotSamplesRelationship(basedir, methylation_data)

### DIFFERENTIAL METHYLATION ANALYSIS ------ ###

print("### Running differential methylation analysis")

# Finding Differentially Methylated Cytosines (DMCs)
if(0 %in% samples_info$Condition & 1 %in% samples_info$Condition) {
  
  if(parameters[7] == 'methylkit') {
    
    # Testing for differential methylation
    dmc_analysis_db <- calculateDiffMeth(methylation_data, overdispersion = "MN", test = "Chisq", mc.cores = parameters[3])
    
    # Converting the methylDiffDB class to a methylDiff object
    dmc_analysis <- as(dmc_analysis_db, "methylDiff")
    
  } else {
    
    # Converting methylation_data to dataframe
    methylation_data <- methylKit::getData(methylation_data)
    
    # log1p transformation
    methylation_data[, 5 : ncol(methylation_data)] <- log1p(methylation_data[, 5 : ncol(methylation_data)])
    
    # Proportional fitting
    library_sizes <- colSums(methylation_data[, seq(5, ncol(methylation_data), 3)])
    methylation_data[, seq(6, ncol(methylation_data), 3)] <- methylation_data[, seq(6, ncol(methylation_data), 3)] * (median(library_sizes) / library_sizes)
    
    # Extract Cs counts
    methylation_data <- methylation_data[, c(1:4, seq(6, ncol(methylation_data), 3))]
    
    # Fold change
    x_cols <- (samples_info$Condition == 0)
    y_cols <- (samples_info$Condition == 1)

    fold_change <- log2(rowMeans(methylation_data[, c(rep(FALSE, 4), x_cols)]) /
                        rowMeans(methylation_data[, c(rep(FALSE, 4), y_cols)]))
    
    # UTest
    pvals <- apply(methylation_data[, 5 : ncol(methylation_data)], 1, function(m) { wilcox.test(m[x_cols], m[y_cols], alternative="two.sided", paired=as.logical(parameters[8]))$p.value })
    
    # FDR
    fdr <- p.adjust(pvals, method="BH")
    
    # Create dataframe
    dmc_analysis <- cbind(methylation_data[, 1:4], pvals, fdr, fold_change)
    colnames(dmc_analysis) <- c("chr", "start", "end", "strand", "pvalue", "qvalue", "meth.diff")
    
    # Converting the dataframe class to a methylDiff object
    dmc_analysis <- as(dmc_analysis, "methylDiff")
    dmc_analysis@sample.ids <- samples_info$SampleID
    dmc_analysis@destranded <- FALSE
    dmc_analysis@assembly <- 'hg18'
    dmc_analysis@context <- 'CpG'
    dmc_analysis@treatment <- samples_info$Condition
    dmc_analysis@resolution <- 'region'
    
  }
  
  # Save to RDS file
  saveRDS(dmc_analysis, paste(basedir, "/DifferentialMethylation_Analysis/", basedir, "_methylDiff.rds", sep=""))
  
  # Deleting database (not needed since the dmc_analysis is saved as an rds file. Also, can be quite large)
  #unlink(paste(basedir, "MethylDB", sep='/'), recursive = TRUE)
  
  # Plot overview of differential methylation events
  plotDMC(basedir, dmc_analysis)
  
  # Extracting differentially methylated regions
  dmc_regions <- extractDMC(basedir, dmc_analysis, parameters[4], parameters[5])
  
  # Plot methylation TSS
  plotMethylationTSS(basedir, dmc_regions, extension=5000)
  
  # Annotating differentially methylated regions
  annotateDMC(basedir, dmc_regions, as.numeric(parameters[9]))
  
}
