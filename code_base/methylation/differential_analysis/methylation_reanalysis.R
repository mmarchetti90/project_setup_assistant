#!/usr/bin/Rscript

# This script analyzes methylation coverage files from Bismark at a base-pair or CpG-island level
# N.B. Current version only works for hg38

### ---------------------------------------- ###

parseArgs <- function() {
  
  # Read command line arguments
  args <- commandArgs()
  
  # Path to rds file with united methylation data
  rds_file <- args[match("--rds_file", args) + 1]
  
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
    
  } else if("--tss_analysis" %in% args) {

    analysis_level <- 'tss'

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
  if("--region_extension" %in% args) {
    
    region_extension <- as.numeric(args[match("--region_extension", args) + 1])
    
  } else {
    
    region_extension <- 0
    
  }
  
  # Islands bp extension upstream and downstream prior to annotation
  if("--correct_batch" %in% args) {
    
    correct_batch <- TRUE
    
  } else {
    
    correct_batch <- FALSE
    
  }

  # DEA file
  if("--dea_file" %in% args) {
    
    dea_file <- args[match("--dea_file", args) + 1]
    
  } else {
    
    dea_file <- ""
    
  }
  
  return(c(rds_file, samples_file, cpus, significance_threshold, significance_parameter, analysis_level, analysis_type, paired, region_extension, correct_batch, dea_file))
  
}

### ---------------------------------------- ###

structureWorkDir <- function(samples_file, analysis_level) {
  
  # Defining basedir name
  if(samples_file == "") {

    if(analysis_level == "island") {

      basedir <- "CpG_Islands_Analysis"

    } else if(analysis_level == "tss") {

      basedir <- "TSS_Analysis"

    } else {

      basedir <- "DMC_Analysis"

    }
    
  } else {

    if(analysis_level == "island") {

      suffix <- "_CpG_Islands_Analysis"

    } else if(analysis_level == "tss") {

      suffix <- "_TSS_Analysis"

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

getSamplesInfo <- function(samples_file) {
  
  info <- read.delim(samples_file, sep = "\t", stringsAsFactors = F)
  info$Condition <- as.integer(as.factor(info$Condition)) - 1
  
  return(info)
  
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
  # png(filename = paste(basedir, "/DifferentialMethylation_Plots/DifferentialMethylationFrequency.png", sep=""), width = 1024, height = 512)
  # diffMethPerChr(dmc_analysis)
  # dev.off()
  
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

annotateDMC_v1 <- function(basedir, dmc_regions, island_extension) {
  
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

### ---------------------------------------- ###

annotateDMC_v2 <- function(basedir, dmc_regions, tss) {
  
  ### N.B. This commented section SHOULD work well, but for some reason the regionCounts function messes up a bit the tss coordinates
  ### so I used the nearest tss
  # Convert dmc_regions and tss to data.frame
  # dmc_regions <- methylKit::getData(dmc_regions)
  # tss <- as.data.frame(tss)
  #
  # # Get coordinates from tss and dmc_regions
  # tss_coords <- paste(tss$seqnames, tss$start, tss$end, sep='_')
  # dmc_coords <- paste(dmc_regions$chr, dmc_regions$start, dmc_regions$end, sep='_')
  # 
  # # Get ordered gene ensembl ids
  # ensembl_ids <- tss$gene_id[match(dmc_coords, tss_coords)]

  # For each region, find closest tss
  ensembl_ids <- tss$gene_id[GenomicRanges::nearest(as(dmc_regions, "GRanges"), tss, ignore.strand=TRUE)]
  
  # Add as first column
  dmc_regions <- methylKit::getData(dmc_regions)
  dmc_regions_cols_order <- c('gene_id', colnames(dmc_regions))
  dmc_regions['gene_id'] <- ensembl_ids
  dmc_regions <- dmc_regions[, dmc_regions_cols_order]
  
  # Saving to file
  write.table(data.frame(dmc_regions), paste(basedir, "/DifferentialMethylation_Analysis/DMC_Annotation.tsv", sep=""), row.names=F, sep="\t", quote=FALSE)
  
}

### ------------------MAIN------------------ ###

library(methylKit)
library(graphics)
library(annotatr)
library(AnnotationHub)
library(org.Hs.eg.db)
library(TxDb.Hsapiens.UCSC.hg38.knownGene)
library(GenomicFeatures)
library(EnrichedHeatmap)
library(ggplot2)
library(circlize)
library(stringr)
library(sva)
library(matrixStats)
library(GenomicRanges)

### LOAD INFO AND STRUCTURE RESULTS DIR ---- ###

print("### Reading info")

# Reading parameters and listing coverage files
parameters <- parseArgs()
print(parameters)

# Structuring workdir
basedir <- structureWorkDir(parameters[2], parameters[6])

### LOAD DATA ----------------------------- ###

# Import samples info
print("### Import samples' info")
samples_info <- getSamplesInfo(parameters[2])

# Import methylation data
print("### Import methylation data")
methylation_data <- readRDS(parameters[1])

# Correction for batch and gender effects
for(var in c("Batch", "Gender")) {
  
  methylation_data <- batch_correction(samples_info, methylation_data, basedir, var, parameters[10])
  
}

# Plotting samples correlation, clustering, and PCA 1,2
print("### Plotting sample relationships")
plotSamplesRelationship(basedir, methylation_data)

# Load DEA
if(parameters[6] == "tss") {

  # Load DEA data, then filter
  print("### Reading DEA file")
  dea <- read.delim(parameters[11], sep='\t', header=TRUE)
  dea <- subset(dea, dea$padj < as.numeric(parameters[4]))

  # Extracting DEGs' tss as Granges objects
  print("### Creating GRange object of DEGs' TSSs")
  dea_ensemble_to_entrez <- mapIds(org.Hs.eg.db, keys=rownames(dea), keytype="ENSEMBL", column = "ENTREZID")
  genes <- genes(TxDb.Hsapiens.UCSC.hg38.knownGene)
  genes <- subset(genes, genes$gene_id %in% as.vector(dea_ensemble_to_entrez))
  tss <- promoters(genes, upstream = as.numeric(parameters[9]), downstream = as.numeric(parameters[9]))

  # Convert TSS entrez ids to ensembl
  tss <- as.data.frame(tss)
  tss$gene_id <- names(dea_ensemble_to_entrez)[match(tss$gene_id, as.vector(dea_ensemble_to_entrez))]
  tss <- as(tss, "GRanges")

  # Converting Granges object chromosome names to match the methylRawList
  tss <- as.data.frame(tss)
  tss$seqnames <- str_replace(tss$seqnames, 'chr', '')
  tss <- as(tss, "GRanges")
  
  # Remove duplicate GRanges
  tss <- as.data.frame(tss)
  tss_coords <- paste(tss$seqnames, tss$start, tss$end, sep='_')
  tss <- tss[! duplicated(tss_coords),]
  tss <- as(tss, "GRanges")

}

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
  
  # Plot overview of differential methylation events
  plotDMC(basedir, dmc_analysis)
  
  # Extracting differentially methylated regions
  dmc_regions <- extractDMC(basedir, dmc_analysis, parameters[4], parameters[5])
  
  if(parameters[6] == "tss") {

    # Annotating differentially methylated regions
    annotateDMC_v2(basedir, dmc_regions, tss)

  } else {

    # Plot methylation TSS
    plotMethylationTSS(basedir, dmc_regions, extension=5000)
  
    # Annotating differentially methylated regions
    annotateDMC_v1(basedir, dmc_regions, as.numeric(parameters[9]))

  }
  
  
}
