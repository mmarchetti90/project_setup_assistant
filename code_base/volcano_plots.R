#!/usr/bin/Rscript

# Volcano plot with highlight of top DEGs

### ---------------------------------------- ###

parseArgs <- function() {
  
  # Read command line arguments
  args <- commandArgs()
  
  # Differential expression file
  deg_file <- args[match("--deg", args) + 1]
  
  # Number of top genes to show
  top_n <- args[match("--top_n", args) + 1]
  
  # Pick top 20 from both upreg and downreg (TRUE) or pick 20 from upreg and 20 from downreg (FALSE, default)?
  if("--absolute_top" %in% args) {
    
    top_toggle <- TRUE
    
  } else {
    
    top_toggle <- FALSE
    
  }
  
  # Adjusted p value threshold
  if("--p_thr" %in% args) {
    
    p_thr <- args[match("--p_thr", args) + 1]
    
  } else {
    
    p_thr <- 0.05
    
  }
  
  # Log2FC threshold
  if("--log2fc_thr" %in% args) {
    
    log2fc_thr <- args[match("--log2fc_thr", args) + 1]
    
  } else {
    
    log2fc_thr <- 0
    
  }
  
  return(c(deg_file, top_n, top_toggle, p_thr, log2fc_thr))
  
}

### ---------------------------------------- ###


### ------------------MAIN------------------ ###

library(ggplot2)
library(RColorBrewer)

### Parse args

parameters <- parseArgs()

### Load DEGs

deg <- read.delim(parameters[1], header = TRUE, sep = "\t", check.names = FALSE)

### Count genes by behaviour

genes_up <- sum((deg$log2FoldChange > as.double(parameters[5])) & (deg$padj < as.double(parameters[4])), na.rm = TRUE)
genes_down <- sum((deg$log2FoldChange < as.double(parameters[5])) & (deg$padj < as.double(parameters[4])), na.rm = TRUE)
genes_ns <- nrow(deg) - genes_up - genes_down

### Prepare data for Volcano plot

# Add log_pval column
deg$log_pval <- - log10(deg$padj)

# Remove genes with NA padj
deg <- subset(deg, ! is.na(deg$padj))

# Prepare color palette
color <- rep(paste("NS (", genes_ns, ")", sep = ""), nrow(deg))
color[(deg$log2FoldChange < - as.double(parameters[5])) & deg$log_pval > - log10(as.double(parameters[4]))] <- paste("Down (", genes_down, ")", sep = "")
color[(deg$log2FoldChange > as.double(parameters[5])) & deg$log_pval > - log10(as.double(parameters[4]))] <- paste("Up (", genes_up, ")", sep = "")
color_palette <- c("green", "red", "gray")
names(color_palette) <- c(paste("Up (", genes_up, ")", sep = ""),
                          paste("Down (", genes_down, ")", sep = ""),
                          paste("NS (", genes_ns, ")", sep = ""))
deg$gene_color <- factor(color)

### Get top genes

if(as.numeric(parameters[2]) > 0) {

  if(parameters[3]) {
    
    deg_sub <- subset(deg, (abs(deg$log2FoldChange) > as.double(parameters[5])) & (deg$padj < as.double(parameters[4])))
    
    genes_to_plot <- deg_sub[order(deg_sub$padj)[1 : min(nrow(deg_sub), as.numeric(parameters[2]))], "gene_symbol"]
    
  } else {
    
    deg_sub_up <- subset(deg, (deg$log2FoldChange > as.double(parameters[5])) & (deg$padj < as.double(parameters[4])))
    
    genes_to_plot <- deg_sub_up[order(deg_sub_up$padj)[1 : min(nrow(deg_sub_up), as.numeric(parameters[2]))], "gene_symbol"]
    
    deg_sub_down <- subset(deg, (deg$log2FoldChange < as.double(parameters[5])) & (deg$padj < as.double(parameters[4])))
    
    genes_to_plot <- c(genes_to_plot, deg_sub_down[order(deg_sub_down$padj)[1 : min(nrow(deg_sub_down), as.numeric(parameters[2]))], "gene_symbol"])
    
  }

} else {

  genes_to_plot <- c()

}

deg$label <- deg$gene_symbol

deg[! deg$gene_symbol %in% genes_to_plot, "label"] <- NA

### Volcano plot

output_name <- "VolcanoPlot.png"

# Main plot
volcano <- ggplot(deg, aes(x = log2FoldChange, y = log_pval, fill = gene_color, label = label)) +
  geom_point(shape = 21, color = "black", size = 3, stroke = 0.5) +
  geom_text(hjust = 0, nudge_x = 0.05) +
  xlab("Log2 Fold Change") +
  ylab("- Log10 Pvalue") +
  scale_fill_manual(values = color_palette) +
  theme(legend.title = element_blank(),
        legend.text = element_text(size=14, face = "bold"),
        axis.text.y = element_text(size=12),
        axis.text.x = element_text(size=12),
        axis.title.y = element_text(size=14, face = "bold"),
        axis.title.x = element_text(size=14, face = "bold"),
        panel.background = element_blank(),
        panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank(),
        panel.border = element_rect(colour = "black", fill=NA, linewidth=3))

# Save to file
ggsave(filename = output_name, plot = volcano, width = 20, height = 25, units = "cm", dpi=300)
