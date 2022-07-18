#BiocManager::install("biocLite")
#source("https://bioconductor.org/biocLite.R")
library(DESeq2)
library(ggplot2)
library(clusterProfiler)
library(biomaRt)
library(ReactomePA)
library(DOSE)
#library(KEGG.db)
library(org.Mm.eg.db)
library(org.Hs.eg.db)

library(pheatmap)
library(genefilter)
library(RColorBrewer)
library(GO.db)
library(topGO)
library(dplyr)
library(gage)
library(ggsci)

# - - - - - - - - - - - - - 
# Import gene counts table
# - - - - - - - - - - - - - 

# Import gene counts table
# - skip first row (general command info)
# - make row names the gene identifiers
countdata <- read.table("./final_counts.txt", header = TRUE, skip = 1, row.names = 1)
#head(countdata)

# Remove .bam + '..' from column identifiers
colnames(countdata) <- gsub(".bam", "", colnames(countdata), fixed = T)
colnames(countdata) <- gsub(".bam", "", colnames(countdata), fixed = T)
colnames(countdata) <- gsub("..", "", colnames(countdata), fixed = T)
colnames(countdata) <- gsub("Aligned.sortedByCoord.out", "", colnames(countdata), fixed = T)

colnames(countdata)
# Remove length/char columns
countdata <- countdata[ ,c(-1:-5)]

# Make sure ID's are correct
head(countdata)


# - - - - - - - - - - - - - 
# Import metadata file
# - - - - - - - - - - - - - 

# Import metadata file
# - make row names the matching sampleID's from the countdata
metadata <- read.delim("./metadata.txt", row.names = 1)

# Add sampleID's to the mapping file
metadata$sampleid <- row.names(metadata)
#head(metadata)
# Reorder sampleID's to match featureCounts column order. 
metadata <- metadata[match(colnames(countdata), metadata$sampleid), ]

# Make sure ID's are correct
head(metadata)

# - - - - - - - - - - - - - 
# Make Deseq2 object
# - - - - - - - - - - - - - 

# - countData : count dataframe
# - colData : sample metadata in the dataframe with row names as sampleID's
# - design : The design of the comparisons to use. 
#            Use (~) before the name of the column variable to compare
ddsMat <- DESeqDataSetFromMatrix(countData = countdata,
                                 colData = metadata,
                                 design = ~Group)

metadata
# Find differential expressed genes
ddsMat <- DESeq(ddsMat)


# Get results from testing with FDR adjust pvalues
results <- results(ddsMat, pAdjustMethod = "fdr", alpha = 0.05)

# Generate summary of testing. 
summary(results)

# Check directionality of the log2 fold changes
## Log2 fold change is set as (LoGlu / HiGlu)
## Postive fold changes = Increased in LoGlu
## Negative fold changes = Decreased in LoGlu
mcols(results, use.names = T)


# - - - - - - - - - - - - - 
# Annotation and write to files
# - - - - - - - - - - - - - 

# Add ENSEMBL
results$ensembl <- mapIds(x = org.Hs.eg.db,
                          keys = row.names(results),
                          column = "ENSEMBL",
                          keytype = "SYMBOL",
                          multiVals = "first")

#make file to save it
new_file <-  as.data.frame(counts(ddsMat))
head(new_file)
new_file$ensembl <- mapIds(x = org.Hs.eg.db,
                           keys = row.names(new_file),
                           column = "ENSEMBL",
                           keytype = "SYMBOL",
                           multiVals = "first")
#new_file
# Write new file to a .txt file
write.table(x = as.data.frame(new_file), 
            file = "results_new_file.txt", 
            sep = '\t', 
            quote = F,
            col.names = NA)


#find out all values where ensembl id is not known
#new_file_final <- new_file[is.na(new_file$ensembl),]
#new_file_final

# Write new_file_final to a .txt file
#write.table(x = as.data.frame(new_file_final), file = "results_new_file_final.txt", sep = '\t', quote = F,col.names = NA)


#THIS IS THE FINAL FILE USED IN BEAVR
install.packages("tidyr")                        # Install & load tidyr package
library("tidyr")
abc <- new_file %>% drop_na()                      # Apply drop_na function
abc

#remove
exportBEAVR <- distinct(abc)

# Write export BEAVR to a .txt file
write.table(x = as.data.frame(exportBEAVR), 
            file = "exportBEAVR.txt", 
            sep = '\t', 
            quote = F,
            col.names = NA)


# Subset for only significant genes (q < 0.05)
results_sig <- subset(results, padj < 0.05)
head(results_sig)
row.names(results_sig)

# Write the annotated results table to a .txt file
write.table(x = as.data.frame(results), 
            file = "results_gene_annotated.txt", 
            sep = '\t', 
            quote = F,
            col.names = NA)

# Write significant annotated results table to a .txt file
write.table(x = as.data.frame(results_sig), 
            file = "results_gene_annotated_significant.txt", 
            sep = '\t', 
            quote = F,
            col.names = NA)

# Write normalized gene counts to a .txt file
write.table(x = as.data.frame(counts(ddsMat), normalized = T), 
            file = 'normalized_counts.txt', 
            sep = '\t', 
            quote = F,
            col.names = NA)

# Write significant normalized gene counts to a .txt file
write.table(x = counts(ddsMat[row.names(results_sig)], normalized = T), 
            file = 'normalized_counts_significant.txt', 
            sep = '\t', 
            quote = F, 
            col.names = NA)

# - - - - - - - - - - - - - 
# Plots
# - - - - - - - - - - - - - 

#PCA

# Convert all samples to rlog
ddsMat_rlog <- rlog(ddsMat, blind = FALSE)

# Plot PCA by column variable
plotPCA(ddsMat_rlog, intgroup = "Group", ntop = 50) +
  theme_bw() + # remove default ggplot2 theme
  geom_point(size = 3) + # Increase point size
  scale_y_continuous(limits = c(-20, 20)) + # change limits to fix figure dimensions
  ggtitle(label = "Principal Component Analysis (PCA)", 
          subtitle = "Top 50 most variable genes") 

#Heatmap
# Convert all samples to rlog
ddsMat_rlog <- rlog(ddsMat, blind = FALSE)

# Gather 30 significant genes and make matrix
mat <- assay(ddsMat_rlog[row.names(results)])[1:40, ]

# Choose which column variables you want to annotate the columns by.
annotation_col = data.frame(
  Group = factor(colData(ddsMat_rlog)$Group), 
  Replicate = factor(colData(ddsMat_rlog)$Replicate),
  row.names = colData(ddsMat_rlog)$sampleid
)

# Specify colors you want to annotate the columns by.
ann_colors = list(
  Group = c(LoGlu = "lightblue", HiGlu = "darkorange"),
  Replicate = c(Rep1 = "darkred", Rep2 = "forestgreen")
)

# Make Heatmap with pheatmap function.
## See more in documentation for customization
pheatmap(mat = mat, 
         color = colorRampPalette(brewer.pal(9, "YlOrBr"))(255), 
         scale = "row", # Scale genes to Z-score (how many standard deviations)
         annotation_col = annotation_col, # Add multiple annotations to the samples
         annotation_colors = ann_colors,# Change the default colors of the annotations
         fontsize = 6.5, # Make fonts smaller
         cellwidth = 55, # Make the cells wider
         show_colnames = F)


#Volcanoplot

# Gather Log-fold change and FDR-corrected pvalues from DESeq2 results
## - Change pvalues to -log10 (1.3 = 0.05)
data <- data.frame(gene = row.names(results),
                   pval = -log10(results$padj), 
                   lfc = results$log2FoldChange)

# Remove any rows that have NA as an entry
data <- na.omit(data)

# Color the points which are up or down
## If fold-change > 0 and pvalue > 1.3 (Increased significant)
## If fold-change < 0 and pvalue > 1.3 (Decreased significant)
data <- mutate(data, color = case_when(data$lfc > 0 & data$pval > 1.3 ~ "Increased",
                                       data$lfc < 0 & data$pval > 1.3 ~ "Decreased",
                                       data$pval < 1.3 ~ "nonsignificant"))

# Make a basic ggplot2 object with x-y values
vol <- ggplot(data, aes(x = lfc, y = pval, color = color))

# Add ggplot2 layers
vol +   
  ggtitle(label = "Volcano Plot", subtitle = "Colored by fold-change direction") +
  geom_point(size = 2.0, alpha = 0.5, na.rm = T) +
  scale_color_manual(name = "Directionality",
                     values = c(Increased = "#008B00", Decreased = "#CD4F39", nonsignificant = "darkgray")) +
  theme_bw(base_size = 14) + # change overall theme
  theme(legend.position = "right") + # change the legend
  xlab(expression(log[2]("Memory_Control" / "Without_treatment"))) + # Change X-Axis label
  ylab(expression(-log[10]("adjusted p-value"))) + # Change Y-Axis label
  geom_hline(yintercept = 1.3, colour = "darkgrey") + # Add p-adj value cutoff line
  scale_y_continuous(trans = "log1p") # Scale yaxis due to large p-values


#MA plot
plotMA(results, ylim = c(-20, 20))
#Plot dispersion
plotDispEsts(ddsMat)

#top gene plot
# Convert all samples to rlog
ddsMat_rlog <- rlog(ddsMat, blind = FALSE)

# Get gene with highest expression
top_gene <- rownames(results)[which.min(results$log2FoldChange)]

# Plot single gene
plotCounts(dds = ddsMat, 
           gene = top_gene, 
           intgroup = "Group", 
           normalized = T, 
           transform = T)
