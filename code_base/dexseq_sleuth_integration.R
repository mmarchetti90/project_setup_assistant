#!/usr/bin/Rscript

# This script reads in differentially expressed transcripts from Sleuth and differentially utilized exons from DEXSeq,
# merges data, and performs enrichment analyses on the genes with concordant changes between the two analysis

### ---------------------------------------- ###

parseArgs <- function() {
  
  # Read command line arguments
  args <- commandArgs()
  
  # Species
  species <- args[match("--species", args) + 1]
  
  # Genome annotation gtf file location
  genome_annotation <- args[match("--genome_annotation", args) + 1]
  
  # Sleuth DEA file location
  sleuth_analysis <- args[match("--sleuth_analysis", args) + 1]
  
  # DEXSeq DEU file location
  dexseq_analysis <- args[match("--dexseq_analysis", args) + 1]
  
  # Adjusted p value threshold
  if("--p_thr" %in% args) {
    
    pval <- args[match("--p_thr", args) + 1]
    
  } else {
    
    pval <- 0.05
    
  }
  
  # msigdbr categories file
  if("--msigdbr_categories_file" %in% args) {
    
    msigdbr_categories_file <- args[match("--msigdbr_categories_file", args) + 1]
    
  } else {
    
    msigdbr_categories_file <- ""
    
  }
  
  return(c(species, genome_annotation, sleuth_analysis, dexseq_analysis, pval, msigdbr_categories_file))
  
}

### ---------------------------------------- ###

loadOrgDB <- function(org) {
  
  if(org == "Homo_sapiens") {
    library(org.Hs.eg.db)
    return(org.Hs.eg.db)
  } else if(org == "Mus_musculus") {
    library(org.Mm.eg.db)
    return(org.Mm.eg.db)
  } else if(org == "Rattus_norvegicus") {
    library(org.Rn.eg.db)
    return(org.Rn.eg.db)
  } else if(org == "Danio_rerio") {
    library(org.Dr.eg.db)
    return(org.Dr.eg.db)
  } else if(org == "Drosophila_melanogaster") {
    library(org.Dm.eg.db)
    return(org.Dm.eg.db)
  } else if(org == "Saccharomyces_cerevisiae") {
    library(org.Sc.eg.db)
    return(org.Sc.eg.db)
  } else if(org == "Caenorhabditis_elegans") {
    library(org.Ce.eg.db)
    return(org.Ce.eg.db)
  } else  {
    print("Unrecognized species...")
    quit(save = "no")
  }
  
}

### ---------------------------------------- ###

loadGenomeAnnotation <- function(gtf_file) {
  
  # Load GTF file
  gtf <- read.table(gtf_file, sep='\t')
  
  # Create transcript_id to gene_id conversion table
  element_info <- gtf[gtf[,3] == "transcript", ncol(gtf)]
  gene_id <- unlist(lapply(element_info, FUN = function(i) { extractElementID(i, "gene_id ") }))
  element_id <- unlist(lapply(element_info, FUN = function(i) { paste(extractElementID(i, "transcript_id "), extractElementID(i, "transcript_version "), sep = '.') }))
  transcript_table <- data.frame("TranscriptID" = element_id,
                                 "GeneID" = gene_id)
  
  # Create exon_id to gene_id conversion table
  element_info <- gtf[gtf[,3] == "exon", ncol(gtf)]
  gene_id <- unlist(lapply(element_info, FUN = function(i) { extractElementID(i, "gene_id ") }))
  element_id <- unlist(lapply(element_info, FUN = function(i) { extractElementID(i, "exon_id ") }))
  exons_table <- data.frame("ExonID" = element_id,
                            "GeneID" = gene_id)
  
  return(list(transcripts = transcript_table, exons = exons_table))
  
}

### ---------------------------------------- ###

extractElementID <- function(string, pattern) {
  
  start_position <- regexpr(pattern, string)[1] + nchar(pattern)
  stop_position <- regexpr(";", substr(string, start_position, nchar(string)))[1] + start_position - 2
  
  element_id <- substr(string, start_position, stop_position)
  
  element_id <- gsub('"', "", element_id)
  
  return(element_id)
  
}

### ---------------------------------------- ###

interpolateData <- function(conv_tbl, trsc, exn) {
  
  # Extracting list of genes whose transcripts are differentially expressed
  transcripts_to_genes <- unique(conv_tbl$transcripts[conv_tbl$transcripts$TranscriptID %in% trsc$target_id, "GeneID"])
  
  # Extracting list of genes whose transcripts are differentially expressed
  exons_to_genes <- unique(conv_tbl$exons[conv_tbl$exons$ExonID %in% exn$featureID, "GeneID"])
  
  # Filtering for common genes
  common_genes <- transcripts_to_genes[transcripts_to_genes %in% exons_to_genes]
  
  return(common_genes)
  
}

### ---------------------------------------- ###

exportGeneList <- function(gtf_file, gl, out) {
  
  # Load GTF file
  gtf <- read.table(gtf_file, sep='\t')
  
  # Create transcript_id to gene_id conversion table
  element_info <- gtf[gtf[,3] == "gene", ncol(gtf)]
  gene_id <- unlist(lapply(element_info, FUN = function(i) { extractElementID(i, "gene_id ") }))
  gene_name <- unlist(lapply(element_info, FUN = function(i) { extractElementID(i, "gene_name ") }))
  annotated_gl <- data.frame("GeneID" = gene_id,
                             "GeneName" = gene_name)
  
  # Selecting only genes in gene_list (gl)
  annotated_gl <- annotated_gl[annotated_gl$GeneID %in% gl,]
  
  # Removing duplicates
  annotated_gl <- unique(annotated_gl)
  
  # Export file
  write.table(annotated_gl, out, row.names = FALSE, sep='\t')
  
}

### ---------------------------------------- ###

runEnricher <- function(out, gl, t2g) {
  
  # Run analysis
  enrichment <- enricher(gene = gl, TERM2GENE = t2g, pvalueCutoff = 0.1, pAdjustMethod = "BH", minGSSize = 10, maxGSSize = 1000)
  
  if(! is.null(enrichment)) {
    
    # Export file
    write.table(enrichment, out, row.names = FALSE, sep='\t')
    
  }
  
  return(enrichment)
  
}

### ---------------------------------------- ###

runEnrichGO <- function(out, odb, ont, gl) {
  
  # Run analysis
  enrichment <- enrichGO(gene = gl, ont = ont, OrgDb = odb, keyType = "ENSEMBL", pvalueCutoff = 0.1, pAdjustMethod = "BH", minGSSize = 10, maxGSSize = 1000)
  
  if(! is.null(enrichment)) {
    
    # Simplifying
    enrichment <- clusterProfiler::simplify(enrichment, cutoff = 0.7)
    
    # Export file
    write.table(enrichment, out, row.names = FALSE, sep='\t')
    
  }
  
  return(enrichment)
  
}

### ------------------MAIN------------------ ###

library(msigdbr)
library(clusterProfiler)
library(enrichplot)

### PARSE ARGS ----------------------------- ###
print("Parsing arguments")

parameters <- parseArgs()
species <- paste(strsplit(parameters[1], "_")[[1]], collapse = " ")
analysis_name <- unlist(strsplit(parameters[3], '/'))
analysis_name <- gsub("Sleuth_LRT_", "", analysis_name[length(analysis_name)])
analysis_name <- gsub(".tsv", "", analysis_name)

### LOADING RESOURCES ---------------------- ###

print("Loading data")

# Loading organism database
org_db <- loadOrgDB(parameters[1])

# Loading gene sets (all of them)
if(species %in% msigdbr_species()$species_name) {
  
  all_gene_sets <- msigdbr(species = species)
  
} else {
  
  print("Unrecognized species...")
  quit(save = "no")
  
}

# Filtering gene sets for user-specified ones
if(parameters[6] != "") {
  
  msigdbr_categories <- as.vector(read.table(parameters[6], sep='\t')[,1])
  all_gene_sets <- all_gene_sets[(all_gene_sets$gs_cat %in% msigdbr_categories) |
                                 (all_gene_sets$gs_subcat %in% msigdbr_categories) |
                                 (all_gene_sets$gs_name %in% msigdbr_categories),]
  
}

# Creating term2gene matrix for enricher, but removing GO terms. They'll be processed and simplified separately
term2gene <- all_gene_sets[!(all_gene_sets$gs_subcat %in% c("GO:BP", "GO:MF", "GO:CC")), c("gs_name", "ensembl_gene")]

# Load gtf genome annotation file and create a transcript_id to gene_id and exon_id to gene_id conversion tables
conversion_tables <- loadGenomeAnnotation(parameters[2])

# Load Sleuth DEA data, filtering for significant elements
transcripts_data <- read.table(parameters[3], sep='\t', header=TRUE, check.names = FALSE)
transcripts_data <- transcripts_data[transcripts_data$qval < as.numeric(parameters[5]),]

# Load DEXSeq DEU data, filtering for significant elements
exon_data <- read.table(parameters[4], sep='\t', header=TRUE, check.names = FALSE)
exon_data <- exon_data[exon_data$padj < as.numeric(parameters[5]),]

### INTERPOLATING DATA --------------------- ###

print("Interpolating data")

# Extracting list of genes with differentially expressed transcripts AND differentially utilized exons
gene_list <- interpolateData(conversion_tables, transcripts_data, exon_data)

# Quit if no genes were found
if(length(gene_list) == 0) {
  
  print("No common genes found...")
  quit(save = "no")
  
}

# Export gene list
output_name <- paste("CommonGenes_", analysis_name, ".tsv", sep = "")
exportGeneList(parameters[2], gene_list, output_name)

### ENRICHMENT ANALYSIS -------------------- ###

print("Running enrichment analyses")

# Init lists of enrichment results to be saved as an RData object
enrichment_analyses <- list()

# Enrichment analysis
output_name <- paste("Enrichment_", analysis_name, ".tsv", sep = "")
enrichment_data <- runEnricher(output_name, gene_list, term2gene)

data_name <- paste("Enrichment_AllSets_", analysis_name, sep = "")
enrichment_analyses[[data_name]] <- enrichment_data

# EnrichGO
for(ontology in c("BP", "MF", "CC")) {
  
  output_name <- paste("Enrichment-GO", ontology, "_", analysis_name, ".tsv", sep = "")
  runEnrichGO(output_name, org_db, ontology, gene_list)
  
  data_name <- paste("Enrichment-GO", ontology, "_", analysis_name, sep = "")
  enrichment_analyses[[data_name]] <- enrichment_data
  
}

# Saving datasets to RData object
output_name <- paste("clusterProfiler_", analysis_name, '.RData.gz', sep = "")
save(parameters, species, all_gene_sets, term2gene, conversion_tables, transcripts_data, exon_data, gene_list, enrichment_analyses, file = output_name, compress = T)
