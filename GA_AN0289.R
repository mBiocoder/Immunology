#=================================================================================================================#
# Exploratory analysis of BD Rhapsody targeted scRNA-Seq from Th subsets (CD25-, CD45RA-, CCR6, CXCR3, CCR4) 
# isolated from healthy skin and blood
# Date: 10.12.2021
# Adapted from https://satijalab.org/seurat/pbmc3k_tutorial.html
# Adapted from https://www.biorxiv.org/content/biorxiv/suppl/2018/03/19/284604.DC1/284604-2.pdf
# Rscript purpose:
# - analyze BD Rhapsody data processed with Seven Bridges Genomics platform
# - dimensionality reduction, clustering and biological annotation of clusters
# - differential expression of each cluster
#=================================================================================================================#
#=================================================================================================================#
# Prepare R session
#=================================================================================================================#

# Load libraries
library(openxlsx)
library(Seurat)
library(plyr)
library(dplyr)
library(SingleR)
library(ggplot2)
library(gridExtra)
library(ComplexHeatmap)

# Define custom theme for ggplot
#theme_set(theme_grey()) # cowplot (loaded along with Seurat v2.3.4 only) sets theme_set(theme_cowplot()) as default
theme_custom <- theme(axis.text.x = element_text(size = 12.8, color = 'black'),
                      axis.text.y = element_text(size = 12.8, color = 'black'), 
                      axis.title.x = element_text(size = 16, margin = margin(6,0,0,0)), 
                      axis.title.y = element_text(size = 16, margin = margin(0,8,0,0)),
                      plot.title = element_text(size = 16, face = 'bold', hjust = 0.5, margin = margin(5,0,10,0)),
                      plot.subtitle = element_text(size = 12, hjust = 0.5, margin = margin(0,0,10,0)), 
                      axis.ticks.x = element_line(size = 0.4, colour = 'black'), 
                      axis.ticks.y = element_line(size = 0.4, colour = 'black'),
                      axis.line = element_blank(), #axis.line = element_line(size = 0.4, colour='black'),
                      panel.grid.minor.y = element_blank(), panel.grid.major.y = element_blank(),
                      panel.grid.minor.x = element_blank(), panel.grid.major.x = element_blank(),
                      panel.border = element_rect(colour = 'black', size = 0.8, fill = NA),
                      panel.background = element_rect(colour = NA, fill = NA),
                      strip.text = element_text(size = 16, face = 'bold', margin = margin(5,0,5,0)), 
                      strip.background = element_blank(), 
                      legend.key = element_blank(), legend.key.size = unit(0.5,'cm'),
                      legend.text = element_text(size = 12.8), 
                      legend.title = element_text(size = 12.8),
                      plot.background = element_rect(colour = NA, fill = NA), aspect.ratio = 1)

#-----------------------------------------------------------------------------------------------------------------#
# !!!DEFINE ALL EXPERIMENT-SPECFIC VARIABLES!!!
#-----------------------------------------------------------------------------------------------------------------#

# Define Experiment_ID according to ngs_sample_list
experimentid1 <- c('EX0011')

# Define Analysis_ID according to ngs_analysis_list
analysisid <- sub('.R', '', stringr::str_split(rstudioapi::getSourceEditorContext()$path, '/', simplify = T)[,5])

# Define |log2fc| cutoff
lfc.cutoff <- 0

# Define analysis description
an.desc <- c('sc_bd_blood_skin_th')
an.descs <- c('sc_th')

#-----------------------------------------------------------------------------------------------------------------#
# Load samples and analysis data
#-----------------------------------------------------------------------------------------------------------------#

# Load table containing the ngs sample list filtered for the data analyzed here
nsl <- read.xlsx('./ngs_sample_list.xlsx', sheet = 1, startRow = 1, colNames = T,
                 rowNames = F, skipEmptyRows = T, skipEmptyCols = T, rows = NULL, cols = NULL, 
                 check.names = T, namedRegion = NULL, na.strings = 'NA', fillMergedCells = F)
nsl <- nsl[nsl$Experiment_ID == experimentid1 ,]#| nsl$Experiment_ID == experimentid2, ]
nsl <- plyr::ddply(nsl, colnames(nsl)[1:18], summarise, sum = NA)

# Update table containing the analyses list to include this analysis
# Include info like which read alignment algorithm was used
# Load, update avoiding duplicates (when script is run several times) and save back overwriting old one
nal <- read.xlsx('./ngs_analysis_list.xlsx', sheet = 1, startRow = 1, colNames = T,
                 rowNames = F, skipEmptyRows = T, skipEmptyCols = T, rows = NULL, cols = NULL, 
                 check.names = T, namedRegion = NULL, na.strings = 'NA', fillMergedCells = F)
nal <- nal[nal$Analysis_ID == analysisid & !is.na(nal$Analysis_ID), ]

# Sample data
metadata <- data.frame(row.names = unique(nsl$Sample_name), 
                       reshape2::colsplit(unique(nsl$Sample_name), '_', c('exp', 'celltype', 'source')))

# Create new folder to store figures from this analysis
dir.create(paste0('figures/', analysisid, '_', an.desc))
dir.create(paste0('figures/', analysisid, '_', an.desc, '/qc'))
#dir.create(paste0('figures/', analysisid, '_', an.desc, '/pca'))
dir.figs <- paste0('figures/', analysisid, '_', an.desc, '/', analysisid, '_', an.descs)
dir.figs.qc <- paste0('figures/', analysisid, '_', an.desc, '/qc/', analysisid, '_', an.descs)
#dir.figs.pca <- paste0('figures/', analysisid, '_', an.desc, '/pca/', analysisid, '_', an.descs)

#=================================================================================================================#
# Analysis of processed data for blood sample
#=================================================================================================================#
#-----------------------------------------------------------------------------------------------------------------#
# Clustering and differential expression analysis with Seurat
#-----------------------------------------------------------------------------------------------------------------#

# Load sample tags
tags <- read.csv(file = './CZ0061_Sample_Tag_Calls.csv', 
                 header = T, skip = 7, row.names = 1)

# Load UMI counts from BD Rhapsody data
bd <- read.csv(file = "./Combined_CZ0061_DBEC_MolsPerCell.csv", 
               header = T, row.names = 1, skip = 7)
colnames(bd) <- sub('HLA.', 'HLA-', colnames(bd))
colnames(bd) <- paste0(stringr::str_split(colnames(bd), '\\.', n = 4, simplify = T)[,1], '_',
                       stringr::str_split(colnames(bd), '\\.', n = 4, simplify = T)[,4])
colnames(bd) <- sub('_Reference_end','_Ref',colnames(bd))
colnames(bd) <- sub('PolyA_','PolyA',colnames(bd))
colnames(bd)[duplicated(colnames(bd))] <- paste0(colnames(bd)[duplicated(colnames(bd))], '2')
colnames(bd) <- sub('_Ref$','',colnames(bd))
colnames(bd)[!(stringr::str_split(colnames(bd), '_', n = 2, simplify = T)[, 1] %in% colnames(bd))] <- 
  sub('_PolyA1','',colnames(bd))[!(stringr::str_split(colnames(bd), '_', n = 2, simplify = T)[, 1] %in% colnames(bd))]
colnames(bd)[!(stringr::str_split(colnames(bd), '_', n = 2, simplify = T)[, 1] %in% colnames(bd))] <- 
  sub('_refseq','',colnames(bd))[!(stringr::str_split(colnames(bd), '_', n = 2, simplify = T)[, 1] %in% colnames(bd))]
#colnames(bd)[!duplicated(stringr::str_split(colnames(bd), '_', n = 2, simplify = T)[, 1])] <- 
#  stringr::str_split(colnames(bd), '_', n = 2, simplify = T)[, 1][
#    !duplicated(stringr::str_split(colnames(bd), '_', n = 2, simplify = T)[, 1])]
#colnames(bd) <- stringr::str_split(colnames(bd), '_', n = 2, simplify = T)[, 1]
#bd <- readRDS(paste0(dir.figs, '_seurat.rds'))
length(unique(colnames(bd)))

# Filter for Th cells only
#bd <- bd[rownames(bd) %in% tags[tags$Sample_Name %in% c('blood_th','skin_th'),]$Cell_Index,]

# Initialize the Seurat object with the raw (non-normalized data).
bd <- CreateSeuratObject(counts = t(bd), min.cells = 0, min.features = 0, project = 'th', meta.data = tags)
VlnPlot(bd, features = c('nFeature_RNA', 'nCount_RNA'), ncol = 2, group.by = 'orig.ident')

# Filter Seurat object for Th cells only
#bd <- subset(bd, subset = !(Sample_Name %in% c('Multiplet','Undetermined')) )
bd <- subset(bd, subset = Sample_Name %in% c('blood_th1','blood_th2','blood_th17','skin_th1','skin_th2','skin_th17'))

# Visualize feature-feature relationships
FeatureScatter(bd, feature1 = 'nCount_RNA', feature2 = 'nFeature_RNA', cols = 'black', group.by = 'orig.ident')

# Filter out cells according to the unique feature counts
bd <- subset(bd, subset = nFeature_RNA >= 20 & nFeature_RNA <= 80 & nCount_RNA <= 2500)
#set.seed(1); bd <- subset(bd, cells = sample(colnames(bd), size = 1000)) #length(colnames(bd)); head(colnames(bd),3)
vp.bd <- VlnPlot(bd, features = c('nFeature_RNA', 'nCount_RNA'), ncol = 2, group.by = 'orig.ident')

# Employ a global-scaling normalization method that normalizes the feature expression measurements for each cell 
# by the total expression, multiplies this by a scale factor (10000 by default), and log-transforms the result
bd <- NormalizeData(bd, normalization.method = 'LogNormalize', scale.factor = 1e4)

# Detection of variable features across the single cells
# Calculate the average expression and dispersion for each feature and z-score for dispersion within each bin
#bd <- FindVariableFeatures(bd, selection.method = 'mean.var.plot', mean.cutoff = c(0.0125, 3), 
#                           dispersion.cutoff = c(0.5, Inf))
bd <- FindVariableFeatures(bd, selection.method = 'vst', nfeatures = length(rownames(bd)), 
                           mean.cutoff = c(0.0125, 3), dispersion.cutoff = c(0.5, Inf))
length(VariableFeatures(bd))
vf.bd <- VariableFeaturePlot(bd, cols = c('black', 'red'), pt.size = 1, log = NULL, assay = NULL) + theme_bw() + theme_custom

# Scaling the data and removing unwanted sources of variation (stored in the scale.data slot)
bd <- ScaleData(bd, features = rownames(bd), vars.to.regress = c('nCount_RNA'))

# Perform linear dimensional reduction
bd <- RunPCA(bd, features = VariableFeatures(bd), verbose = F)

# Examine and visualize PCA results a few different ways (VizDimReduction, DimPlot, and DimHeatmap)
print(bd[['pca']], dims = 1:5, nfeatures = 5, projected = F)
VizDimLoadings(bd, dims = c(1,2))
DimPlot(bd, reduction = 'pca')

# Scores each feature in the dataset (including features not included in the PCA) based on their correlation with the PCs
bd <- ProjectDim(bd, do.center = T, verbose = F)
print(bd[['pca']], dims = 1:5, nfeatures = 5, projected = T)
VizDimLoadings(bd, dims = c(1,2), projected = T)
DimPlot(bd, projected = T, reduction = 'pca')

# Explore the primary sources of heterogeneity to decide which PCs to include for further downstream analyses
DimHeatmap(bd, dims = 1, cells = NULL, balanced = T)
#DimHeatmap(bd, dims = 1:20, cells = 500, balanced = T)
hm.pc.bd <- DimHeatmap(bd, dims = 1:30, cells = 500, balanced = T, fast = F)# + ggtitle('Heatmap of PC markers from 500 cells (blood)')
hm.pc.bd2 <- DimHeatmap(bd, dims = 30:50, cells = 500, balanced = T, fast = F)# + ggtitle('Heatmap of PC markers from 500 cells (blood)')

# Determine statistically significant principal components
bd <- JackStraw(bd, dims = 50, prop.freq = 0.1, num.replicate = 100, maxit = 1000)
bd <- ScoreJackStraw(bd, dims = 1:50)
JackStrawPlot(bd, dims = 1:50)
JackStrawPlot(bd, dims = 15:25)
ep.bd <- ElbowPlot(bd, ndims = 50)

# Cluster cells (the ‘granularity’ (number of clusters) of the downstream clustering is set to 0.6-1.2 for 3K cells)
bd <- FindNeighbors(bd, reduction = 'pca', k.param = 30, dims = 1:14, do.plot = F)
bd <- FindClusters(bd, resolution = 1.2, algorithm = 1, n.start = 10, n.iter = 10, random.seed = 0)
#DimPlot(bd, reduction = 'tsne', dims = c(1,2), group.by = 'seurat_clusters', split.by = 'orig.ident', label = T, label.size = 4, pt.size = 1.5, cols = NULL) + theme_bw() + theme_custom + xlab('tSNE 1') + ylab('tSNE 2') + scale_y_continuous(breaks = seq(-25,25,25)) + scale_x_continuous(breaks = seq(-25,25,25)) #+ NoLegend()
#DimPlot(bd, reduction = 'tsne', dims = c(1,2), group.by = 'Sample_Name', split.by = 'orig.ident', label = T, label.size = 4, pt.size = 1.5, cols = NULL) + theme_bw() + theme_custom + xlab('tSNE 1') + ylab('tSNE 2') + scale_y_continuous(breaks = seq(-25,25,25)) + scale_x_continuous(breaks = seq(-25,25,25)) #+ NoLegend()

# Run Non-linear dimensional reduction (tSNE and UMAP)
bd <- RunTSNE(bd, reduction = 'pca', dims = 1:14, seed.use = 1, tsne.method = 'Rtsne', dim.embed = 2, perplexity = 30)
DimPlot(bd, reduction = 'tsne', dims = c(1,2), group.by = 'RNA_snn_res.1.2', split.by = 'orig.ident', label = T, label.size = 8, pt.size = 1.5, cols = NULL) + theme_bw() + theme_custom
#bd <- RunUMAP(bd, reduction = 'pca', dims = 1:14, seed.use = 3, spread = 10, min.dist = 0.001, reduction.name = 'umap')
bd <- RunUMAP(bd, reduction = 'pca', dims = 1:14, seed.use = 3, spread = 5, min.dist = 0.3, reduction.name = 'umap')
DimPlot(bd, reduction = 'umap', dims = c(1,2), group.by = 'RNA_snn_res.1.2', split.by = 'orig.ident', label = T, label.size = 8, pt.size = 1.5, cols = NULL) + theme_bw() + theme_custom
#ts.sp.bd <- DimPlot(bd, reduction = 'tsne', dims = c(1,2), group.by = 'RNA_snn_res.1.2', split.by = 'RNA_snn_res.0.4', 
#                    label = F, label.size = 8, pt.size = 1, cols = NULL) + scale_y_continuous(breaks = c(-25,0,25)) + theme_bw() + theme_custom
vp.bdc <- VlnPlot(bd, features = c('nFeature_RNA', 'nCount_RNA'), ncol = 2, group.by = 'RNA_snn_res.1.2')
#vp.bdc.sum <- VlnPlot(bd, features = c('nFeature_RNA', 'nCount_RNA'), ncol = 2, group.by = 'RNA_snn_res.0.4')

# Finding differentially expressed features (cluster biomarkers)
#Idents(bd) <- bd$RNA_snn_res.1.2
bd.markers <- FindAllMarkers(bd, min.pct = 0.1, logfc.threshold = 0.25, min.diff.pct = -Inf, only.pos = F, 
                             test.use = 'wilcox', return.thresh = 0.01, max.cells.per.ident = Inf, random.seed = 1)
#bd.markers <- bd.markers[order(bd.markers$cluster, decreasing = F),]
write.csv(bd.markers, file = paste0(dir.figs, '_cluster_markers.csv'), row.names = T)
#write.csv(bd.markers.sum, file = paste0(dir.figs, '_cluster_markers_summarized.csv'), row.names = T)
#bd.markers.sum <- read.csv(file = paste0(dir.figs, '_cluster_markers_summarized.csv'), row.names = 1)
bd.top <- bd.markers[bd.markers$avg_log2FC >= 0, ] %>% group_by(cluster) %>% top_n(n = 10, wt = avg_log2FC)
#colnames(bd.markers)
gn <- c('S1PR1', 'SELL', 'CD69', 'LGALS3', 'RUNX3', 'FABP5')

# Visualization of marker genes expression (VlnPlot, FeaturePlot, RidgePlot, CellScatter, and DotPlot)
#VlnPlot(bd, features = bd.top$gene[1:6], cols = NULL, pt.size = 1, adjust = 1,  y.max = NULL, ncol = NULL, slot = 'data')
#VlnPlot(bd, features = gn, cols = NULL, pt.size = 1, adjust = 1,  y.max = NULL, ncol = NULL, slot = 'data')
FeaturePlot(bd, features = bd.top$gene[1], dims = c(1,2), cols = c('lightgrey', 'blue'), pt.size = NULL, 
            min.cutoff = NA, max.cutoff = 5, reduction = 'tsne', blend = FALSE, blend.threshold = 0.5, 
            order = NULL, label = F, label.size = 4, ncol = NULL, coord.fixed = T)
hm.mk.bd <- DoHeatmap(bd, features = bd.top$gene, disp.min = -2.5, disp.max = NULL, slot = 'scale.data', label = T, 
                      size = 5.5, hjust = 0, angle = 0, combine = T, group.by = 'RNA_snn_res.1.2', assay = 'RNA') + 
  ggtitle('Heatmap of top 10 cluster markers (Th cells)')# + NoLegend()
#saveRDS(bd, file = paste0(dir.figs, '_seurat.rds'))
bd <- readRDS('./GA_AN0289_sc_th_seurat.rds')
bd.markers <- read.csv(file = paste0(dir.figs, '_cluster_markers.csv'), row.names = 1)

#=================================================================================================================#
# Cluster annotation
#=================================================================================================================#
#-----------------------------------------------------------------------------------------------------------------#
# Biological annotation of clusters from blood sample using SingleR
#-----------------------------------------------------------------------------------------------------------------#

train_single <- trainSingleR(GetAssayData(bd, slot = 'data'),bd@active.ident, genes = 'de',  sd.thresh = 1,
                      de.method = c("classic", "wilcox", "t"),
                      de.n = NULL,
                      de.args = list(),
                      aggr.ref = FALSE,
                      aggr.args = list(),
                      recompute = TRUE,
                      restrict = NULL,
                      assay.type = "logcounts",
                      check.missing = TRUE  )

sd.bd <- classifySingleR(
  GetAssayData(bd, slot = 'data'),
  train_single,
  quantile = 0.8,
  fine.tune = TRUE,
  tune.thresh = 0.05,
  sd.thresh = NULL,
  prune = TRUE,
  assay.type = "logcounts",
  check.missing = TRUE)

sr.bd <- readRDS('GA_AN0289_sc_th_singler.rds')

# Annotate SingleR object with tSNE coordinates and clusters from Seurat object
#sr.bd$seurat = seurat.object # (optional)
sr.bd$meta.data$orig.ident <- bd@meta.data$orig.ident # the original identities, if not supplied in 'annot'
sr.bd$meta.data$xy <- bd@reductions$tsne@cell.embeddings # the tSNE coordinates (Embeddings(bd, reduction = 'tsne'))
sr.bd$meta.data$xy.um <- bd@reductions$umap@cell.embeddings # the tSNE coordinates (Embeddings(bd, reduction = 'umap'))
sr.bd$meta.data$clusters <- bd$RNA_snn_res.1.2 # the Seurat clusters, if 'clusters' not provided (Idents(bd))
#sr.bd$singler[[2]]$SingleR.single$labels <- plyr::mapvalues(sr.bd$singler[[2]]$SingleR.single$labels, 
#                                                            from = levels(factor(sr.bd$singler[[2]]$SingleR.single$labels)), 
#                                                            to = c('CD4 Tn', 'CD4 Tcm', 'CD4 Tem', 
#                                                                   'CD8 Tn', 'CD8 Tcm', 'CD8 Tem', 
#                                                                   'NK cell', 'Treg'))

## Correlation of variable genes between the experiment and the reference dataset
#SingleR.DrawScatter(sc_data = GetAssayData(bd, slot = 'data'), cell_id = 1, ref = hpca, sample_id = 1)
#
## Correlation between a single cell from the experiment and all reference cells
#SingleR.DrawBoxPlot(sc_data = GetAssayData(bd, slot = 'data'), cell_id = 1, ref = blueprint_encode, 
#                    main_types = T, labels.use = NULL)$p
#
## Transfer labels to Seurat object
#bd[['celltype']] <- sr.bd$singler[[2]]$SingleR.single$labels[match(sr.bd$singler[[2]]$SingleR.single$cell.names, colnames(bd))]
##saveRDS(bd, file = paste0(dir.figs, '_seurat.rds'))
#
## Number of cells per annotated cluster
#sr.bdn1 <- table(sr.bd$meta.data$orig.ident, sr.bd$meta.data$clusters); knitr::kable(sr.bdn1)
#sr.bdn2 <- table(sr.bd$singler[[2]]$SingleR.single$labels, sr.bd$meta.data$orig.ident); knitr::kable(sr.bdn2)
#sr.bdn3 <- table(sr.bd$singler[[2]]$SingleR.single$labels, sr.bd$meta.data$clusters); knitr::kable(sr.bdn3)
#sr.bdn4 <- table(sr.bd$singler[[2]]$SingleR.single.main$labels, sr.bd$meta.data$orig.ident); knitr::kable(sr.bdn4)
#sr.bdn5 <- table(sr.bd$meta.data$orig.ident, sr.bd$meta.data$clusters.sum); knitr::kable(sr.bdn5)
#sr.bdn6 <- table(sr.bd$singler[[2]]$SingleR.single$labels, sr.bd$meta.data$orig.ident); knitr::kable(sr.bdn6)
#sr.bdn7 <- table(sr.bd$singler[[2]]$SingleR.single$labels, sr.bd$meta.data$clusters.sum); knitr::kable(sr.bdn7)
#sr.bdn8 <- table(sr.bd$singler[[2]]$SingleR.single.main$labels, sr.bd$meta.data$orig.ident); knitr::kable(sr.bdn8)
#
## Heatmap of the aggregated scores before fine-tuning for the main cell types:
#SingleR.DrawHeatmap(sr.bd$singler[[2]]$SingleR.single, top.n = Inf, clusters = sr.bd$meta.data$cluster)#orig.ident)
#
## Define colors
#levels(factor(sr.bd$meta.data$clusters))
#levels(factor(sr.bd$singler[[2]]$SingleR.single$labels[,1]))
#col1 <- c('darkorange2', 'firebrick1', 'darkgoldenrod1', 'steelblue1', 'grey70', 'plum3', 'darkorchid2', 
#          'seagreen3', 'grey40', 'royalblue1', 'olivedrab', 'orchid1', 'wheat2')
##col2 <- c('darkgoldenrod1', 'firebrick1', 'darkorange2', 'steelblue1', 'royalblue1', 'plum3', 
##          'darkorchid2', 'olivedrab')
#col2 <- c('firebrick1', 'grey70', 'darkgoldenrod1', 'royalblue1', 'plum3', 'steelblue1', 'darkorchid2', 'olivedrab')
#col3 <- c('darkorange2', 'firebrick1', 'olivedrab', 'grey40', 'firebrick1', 'darkgoldenrod1', 'darkorange2', 
#          'darkorchid2', 'steelblue1', 'royalblue1', 'seagreen3', 'darkorchid2', 'darkorange2')
#
## Draw tSNE plots swith cluster from Seurat object
#ts.bd <- SingleR.PlotTsne(sr.bd$singler[[2]]$SingleR.single, sr.bd$meta.data$xy, do.labels = T, do.letters = F, col = col1, 
#                          labels = sr.bd$meta.data$clusters, label.size = 4, dot.size = 1, title = 'blood')$p + 
#  scale_color_manual(values = col1, guide = guide_legend(override.aes = list(size = 4, alpha = 0.8))) + 
#  scale_y_continuous(breaks = seq(-20,20,20)) + scale_x_continuous(breaks = seq(-25,25,25)) + theme_custom
#ts.bd.sum <- SingleR.PlotTsne(sr.bd$singler[[2]]$SingleR.single, sr.bd$meta.data$xy, do.labels = T, do.letters = F, col = col1, 
#                              labels = sr.bd$meta.data$clusters.sum, label.size = 4, dot.size = 1, title = 'blood')$p + 
#  scale_color_manual(values = col1, guide = guide_legend(override.aes = list(size = 4, alpha = 0.8))) + 
#  scale_y_continuous(breaks = seq(-20,20,20)) + scale_x_continuous(breaks = seq(-25,25,25)) + theme_custom
## Draw tSNE plots with annotation from SingleR at single cell level
#ta.bd <- SingleR.PlotTsne(sr.bd$singler[[2]]$SingleR.single, sr.bd$meta.data$xy, do.labels = T, do.letters = F, col = col2, 
#                          labels = sr.bd$singler[[2]]$SingleR.single$labels, label.size = 4, dot.size = 1, title = 'blood')$p + 
#  scale_color_manual(values = col2, guide = guide_legend(override.aes = list(size = 4, alpha = 0.8))) + 
#  scale_y_continuous(breaks = seq(-20,20,20)) + scale_x_continuous(breaks = seq(-25,25,25)) + theme_custom
#ta.bd.nl <- SingleR.PlotTsne(sr.bd$singler[[2]]$SingleR.single, sr.bd$meta.data$xy, do.labels = F, do.letters = F, col = col2, 
#                             labels = sr.bd$singler[[2]]$SingleR.single$labels, label.size = 4, dot.size = 1, title = 'blood')$p + 
#  scale_color_manual(values = col2, guide = guide_legend(override.aes = list(size = 4, alpha = 0.8))) + 
#  scale_y_continuous(breaks = seq(-20,20,20)) + scale_x_continuous(breaks = seq(-25,25,25)) + theme_custom
## Draw tSNE plots with annotation from SingleR at cluster level
#tc.bd <- SingleR.PlotTsne(sr.bd$singler[[2]]$SingleR.clusters, sr.bd$meta.data$xy, do.labels = T, do.letters = F, col = col3, 
#                          clusters = sr.bd$meta.data$clusters, #singlec$singler[[2]]$SingleR.single$clusters$cl, 
#                          labels = sr.bd$singler[[2]]$SingleR.clusters$labels, label.size = 4, dot.size = 1, title = 'blood')$p + 
#  scale_color_manual(values = col3, guide = guide_legend(override.aes = list(size = 4, alpha = 0.8))) + 
#  scale_y_continuous(breaks = seq(-20,20,20)) + scale_x_continuous(breaks = seq(-25,25,25)) + theme_custom
#tc.bd.nl <- SingleR.PlotTsne(sr.bd$singler[[2]]$SingleR.clusters, sr.bd$meta.data$xy, do.labels = F, do.letters = F, col = col3, 
#                             clusters = sr.bd$meta.data$clusters, #singlec$singler[[2]]$SingleR.single$clusters$cl, 
#                             labels = sr.bd$singler[[2]]$SingleR.clusters$labels, label.size = 4, dot.size = 1, title = 'blood')$p + 
#  scale_color_manual(values = col3, guide = guide_legend(override.aes = list(size = 4, alpha = 0.8))) + 
#  scale_y_continuous(breaks = seq(-20,20,20)) + scale_x_continuous(breaks = seq(-25,25,25)) + theme_custom
## Draw UMAP plots swith cluster from Seurat object
#um.bd <- SingleR.PlotTsne(sr.bd$singler[[2]]$SingleR.single, sr.bd$meta.data$xy.um, do.labels = T, do.letters = F, col = col1, 
#                          labels = sr.bd$meta.data$clusters, label.size = 4, dot.size = 1, title = 'blood')$p + 
#  xlab('UMAP 1') + ylab('UMAP 2') + scale_y_continuous(breaks = seq(-25,25,25)) + scale_x_continuous(breaks = seq(-25,25,25)) + 
#  scale_color_manual(values = col1, guide = guide_legend(override.aes = list(size = 4, alpha = 0.8))) + theme_custom
#um.bd.sum <- SingleR.PlotTsne(sr.bd$singler[[2]]$SingleR.single, sr.bd$meta.data$xy.um, do.labels = T, do.letters = F, col = col1, 
#                              labels = sr.bd$meta.data$clusters.sum, label.size = 4, dot.size = 1, title = 'blood')$p + 
#  xlab('UMAP 1') + ylab('UMAP 2') + scale_y_continuous(breaks = seq(-25,25,25)) + scale_x_continuous(breaks = seq(-25,25,25)) + 
#  scale_color_manual(values = col1, guide = guide_legend(override.aes = list(size = 4, alpha = 0.8))) + theme_custom
## Draw UMAP plots with annotation from SingleR at single cell level
#ua.bd <- SingleR.PlotTsne(sr.bd$singler[[2]]$SingleR.single, sr.bd$meta.data$xy.um, do.labels = T, do.letters = F, col = col2, 
#                          labels = sr.bd$singler[[2]]$SingleR.single$labels, label.size = 4, dot.size = 1, title = 'blood')$p + 
#  xlab('UMAP 1') + ylab('UMAP 2') + scale_y_continuous(breaks = seq(-25,25,25)) + scale_x_continuous(breaks = seq(-25,25,25)) + 
#  scale_color_manual(values = col2, guide = guide_legend(override.aes = list(size = 4, alpha = 0.8))) + theme_custom
#ua.bd.nl <- SingleR.PlotTsne(sr.bd$singler[[2]]$SingleR.single, sr.bd$meta.data$xy.um, do.labels = F, do.letters = F, col = col2, 
#                             labels = sr.bd$singler[[2]]$SingleR.single$labels, label.size = 4, dot.size = 1, title = 'blood')$p + 
#  xlab('UMAP 1') + ylab('UMAP 2') + scale_y_continuous(breaks = seq(-25,25,25)) + scale_x_continuous(breaks = seq(-25,25,25)) + 
#  scale_color_manual(values = col2, guide = guide_legend(override.aes = list(size = 4, alpha = 0.8))) + theme_custom
## Draw UMAP plots with annotation from SingleR at cluster level
#uc.bd <- SingleR.PlotTsne(sr.bd$singler[[2]]$SingleR.clusters, sr.bd$meta.data$xy.um, do.labels = T, do.letters = F, col = col3, 
#                          clusters = sr.bd$meta.data$clusters, #singlec$singler[[2]]$SingleR.single$clusters$cl, 
#                          labels = sr.bd$singler[[2]]$SingleR.clusters$labels, label.size = 4, dot.size = 1, title = 'blood')$p + 
#  xlab('UMAP 1') + ylab('UMAP 2') + scale_y_continuous(breaks = seq(-25,25,25)) + scale_x_continuous(breaks = seq(-25,25,25)) + 
#  scale_color_manual(values = col3, guide = guide_legend(override.aes = list(size = 4, alpha = 0.8))) + theme_custom
#uc.bd.nl <- SingleR.PlotTsne(sr.bd$singler[[2]]$SingleR.clusters, sr.bd$meta.data$xy.um, do.labels = F, do.letters = F, col = col3, 
#                             clusters = sr.bd$meta.data$clusters, #singlec$singler[[2]]$SingleR.single$clusters$cl, 
#                             labels = sr.bd$singler[[2]]$SingleR.clusters$labels, label.size = 4, dot.size = 1, title = 'blood')$p + 
#  xlab('UMAP 1') + ylab('UMAP 2') + scale_y_continuous(breaks = seq(-25,25,25)) + scale_x_continuous(breaks = seq(-25,25,25)) + 
#  scale_color_manual(values = col3, guide = guide_legend(override.aes = list(size = 4, alpha = 0.8))) + theme_custom

#-----------------------------------------------------------------------------------------------------------------#
# Cluster annotation and gene expression visualization
#-----------------------------------------------------------------------------------------------------------------#

# Prepare data table
df.bd <- data.frame(tx = bd@reductions$tsne@cell.embeddings[,1],
                    ty = bd@reductions$tsne@cell.embeddings[,2],
                    ux = bd@reductions$umap@cell.embeddings[,1],
                    uy = bd@reductions$umap@cell.embeddings[,2], 
                    cells = colnames(bd), 
                    tissue = stringr::str_split(bd$Sample_Name, "_", n = 2, simplify = T)[,2], 
                    subset = stringr::str_split(bd$Sample_Name, "_", n = 2, simplify = T)[,1], 
                    tissue_subset = sub('_',' ',bd$Sample_Name), 
                    cluster = bd$RNA_snn_res.1.2, 
                    t(as.matrix(GetAssayData(bd, slot = 'data'))))#[gn, ])))
df.bd <- reshape2::melt(df.bd, id.vars = c('tx','ty', 'ux','uy','cells','tissue','subset','tissue_subset','cluster'))
df.bd$variable <- sub('\\.','-',df.bd$variable)
df.bd$tissue_subset <- factor(df.bd$tissue_subset, levels = levels(factor(df.bd$tissue_subset))[c(1,3,5,2,4,6)])

# Visualization of clusters
col <- list(cluster = c(RColorBrewer::brewer.pal(n = 12, name = 'Paired')[c(1:12)], 'grey60'))
ts.cl <- ggplot(df.bd[!duplicated(df.bd$cells),], aes(x = tx, y = ty, color = cluster, fill = cluster)) + 
  geom_point(shape = 16, size = 1.5, alpha = 0.6) + 
  scale_x_continuous(breaks = seq(-20,20,20)) + scale_y_continuous(breaks = seq(-20,20,20)) + 
  #scale_color_manual(values = col$cluster) + guides(color = guide_legend(override.aes = list(size = 4, alpha = 0.6))) + 
  scale_color_manual(values = col$cluster, guide = guide_legend(override.aes = list(size = 4, alpha = 0.6))) + 
  scale_fill_manual(values = col$cluster) + 
  theme_custom + xlab('tSNE 1') + ylab('tSNE 2') + ggtitle('blood and skin Th cells')#, subtitle = 'Louvain clusters')
um.cl <- ggplot(df.bd[!duplicated(df.bd$cells),], aes(x = ux, y = uy, color = cluster, fill = cluster)) + 
  geom_point(shape = 16, size = 1.5, alpha = 0.6) + 
  scale_x_continuous(breaks = seq(-30,30,15)) + scale_y_continuous(breaks = seq(-15,15,15)) + 
  scale_color_manual(values = col$cluster, guide = guide_legend(override.aes = list(size = 4, alpha = 0.6))) + 
  scale_fill_manual(values = col$cluster) + 
  theme_custom + xlab('UMAP 1') + ylab('UMAP 2') + ggtitle('blood and skin Th cells')

# Visualization of subsets
#col$bw1 <- rep(c('grey70', 'grey40'),3)
#col$bw2 <- rep(c('white', 'black'),3)
col$ts <- RColorBrewer::brewer.pal(n = 12, name = 'Paired')[c(1,3,5,2,4,6)]
ts.ct <- ggplot(df.bd[!duplicated(df.bd$cells),], 
                aes(x = tx, y = ty, color = tissue_subset, fill = tissue_subset)) + 
  geom_point(shape = 21, size = 1.5, alpha = 0.8) + 
  scale_x_continuous(breaks = seq(-20,20,20)) + scale_y_continuous(breaks = seq(-20,20,20)) + 
  scale_color_manual(values = col$ts) + guides(color = guide_legend(override.aes = list(size = 4, alpha = 0.8))) + 
  scale_fill_manual(values = col$ts) + 
  #scale_shape_manual(values = c(21,21,22,22,24,24)) + 
  theme_custom + xlab('tSNE 1') + ylab('tSNE 2') + ggtitle('blood and skin Th cells')
um.ct <- ggplot(df.bd[!duplicated(df.bd$cells),], aes(x = ux, y = uy, color = tissue_subset, fill = tissue_subset)) + 
  geom_point(shape = 21, size = 1.5, alpha = 0.8) + 
  scale_x_continuous(breaks = seq(-30,30,15)) + scale_y_continuous(breaks = seq(-15,15,15)) + 
  scale_color_manual(values = col$ts) + guides(color = guide_legend(override.aes = list(size = 4, alpha = 0.8))) + 
  scale_fill_manual(values = col$ts) + 
  theme_custom + xlab('UMAP 1') + ylab('UMAP 2') + ggtitle('blood and skin Th cells')

# Visualization of marker gene expression
#df.bd[df.bd$variable %in% head(bd.markers[order(bd.markers$p_val_adj, decreasing = F),]$gene,10), ]
#col$hm <- colorRampPalette(c('grey80', 'hotpink3', 'royalblue4'))(50)
#col$hm <- viridis::viridis_pal()(40)
col$hm <- colorRampPalette(c('grey85', 'bisque2', 'hotpink3', 'royalblue4'))(50)
ts.mk.bd <- ggplot(df.bd[df.bd$variable %in% (bd.markers[bd.markers$avg_log2FC >= 0, ] %>% 
                           group_by(cluster) %>% top_n(n = 2, wt = avg_log2FC))$gene, ], 
                   aes(x = tx, y = ty, color = value)) + geom_point(shape = 16, size = 1, alpha = 1) + 
  scale_x_continuous(breaks = seq(-20,20,20)) + scale_y_continuous(breaks = seq(-20,20,20)) + 
  #scale_color_gradient(expression(''*log[2]*'(exp)'), low = 'gray', high = 'darkblue', na.value = 'darkblue', limits = c(0,8)) + 
  #guides(color = guide_colorbar(title.position = 'left') + 
  scale_color_gradientn(expression(''*log[2]*'(exp)'), colours = col$hm, oob = scales::squish, limits = c(0,8)) + 
  facet_wrap(~variable, ncol = 5) + theme_custom + xlab('tSNE 1') + ylab('tSNE 2') + ggtitle('blood and skin Th cells')
um.mk.bd <- ggplot(df.bd[df.bd$variable %in% (bd.markers[bd.markers$avg_log2FC >= 0, ] %>% 
                                                group_by(cluster) %>% top_n(n = 2, wt = avg_log2FC))$gene, ], 
                   aes(x = ux, y = uy, color = value)) + geom_point(shape = 16, size = 1, alpha = 1) + 
  scale_x_continuous(breaks = seq(-30,30,15)) + scale_y_continuous(breaks = seq(-15,15,15)) + 
  scale_color_gradientn(expression(''*log[2]*'(exp)'), colours = col$hm, oob = scales::squish, limits = c(0,8)) + 
  facet_wrap(~variable, ncol = 5) + theme_custom + xlab('UMAP 1') + ylab('UMAP 2') + ggtitle('blood and skin Th cells')

# Visualization of gene expression for all genes
#ts.all.bd <- ggplot(df.bd, aes(x = tx, y = ty, color = value)) + geom_point(shape = 16, size = 1, alpha = 1) + 
#  scale_x_continuous(breaks = seq(-15,15,15)) + scale_y_continuous(breaks = seq(-20,20,20)) + 
#  scale_color_gradientn('', colours = col$hm, oob = scales::squish, limits = c(0,8)) + 
#  facet_wrap(~variable, ncol = 30) + theme_custom + xlab('tSNE 1') + ylab('tSNE 2') + 
#  theme(axis.title.x = element_text(size = 14, margin = margin(6,0,0,0)), 
#        axis.title.y = element_text(size = 14, margin = margin(0,8,0,0)),
#        plot.title = element_text(size = 14, face = 'bold', hjust = 0.5, margin = margin(5,0,10,0)),
#        strip.text = element_text(size = 14, face = 'bold', margin = margin(5,0,5,0)))
ts.all.bd <- Map(x=unique(df.bd$variable), function(x) {
  qt <- ceiling(quantile(df.bd[df.bd$variable == x,]$value, 0.99))
  lm <- ifelse(qt<=3,3, ifelse(qt<=6,6, ifelse(qt<=9,9)))
  ggplot(df.bd[df.bd$variable %in% x,], aes(x = tx, y = ty, color = value)) + 
    geom_point(shape = 16, size = 1.5, alpha = 1) + 
    scale_x_continuous(breaks = seq(-20,20,20)) + scale_y_continuous(breaks = seq(-20,20,20)) + 
    scale_color_gradientn(expression(''*log[2]*'(exp)'), colours = col$hm, oob = scales::squish, 
                          limits = c(0, lm), breaks = seq(0, lm, length.out = 4)) + 
    facet_wrap(~variable, ncol = 30) + theme_custom + xlab('tSNE 1') + ylab('tSNE 2') + 
    ggtitle('blood and skin Th cells') + 
    theme(axis.title.x = element_text(size = 14, margin = margin(6,0,0,0)), 
          axis.title.y = element_text(size = 14, margin = margin(0,8,0,0)),
          plot.title = element_text(size = 14, face = 'bold', hjust = 0.5, margin = margin(5,0,10,0)),
          strip.text = element_text(size = 14, face = 'bold', margin = margin(5,0,5,0))) })
#gridExtra::marrangeGrob(grobs = ts.all.bd, nrow = 2, ncol = 2)
um.all.bd <- Map(x=unique(df.bd$variable), function(x)  {
  qt <- ceiling(quantile(df.bd[df.bd$variable == x,]$value, 0.99))
  lm <- ifelse(qt<=3,3, ifelse(qt<=6,6, ifelse(qt<=9,9)))
  ggplot(df.bd[df.bd$variable %in% x,], aes(x = ux, y = uy, color = value)) + 
    geom_point(shape = 16, size = 1.5, alpha = 1) + 
    scale_x_continuous(breaks = seq(-30,30,15)) + scale_y_continuous(breaks = seq(-15,15,15)) + 
    scale_color_gradientn(expression(''*log[2]*'(exp)'), colours = col$hm, oob = scales::squish, 
                          limits = c(0, lm), breaks = seq(0, lm, length.out = 4)) + 
    facet_wrap(~variable, ncol = 30) + theme_custom + xlab('UMAP 1') + ylab('UMAP 2') + 
    ggtitle('blood and skin Th cells') + 
    theme(axis.title.x = element_text(size = 14, margin = margin(6,0,0,0)), 
          axis.title.y = element_text(size = 14, margin = margin(0,8,0,0)),
          plot.title = element_text(size = 14, face = 'bold', hjust = 0.5, margin = margin(5,0,10,0)),
          strip.text = element_text(size = 14, face = 'bold', margin = margin(5,0,5,0))) })

# Visualization of gene expression: TRC chains
#ts.mk.bd.tcr <- ggplot(df.bd[df.bd$variable %in% rownames(bd)[grep('TR', rownames(bd))][c(1,3,4)],], aes(x = tx, y = ty, color = value)) + geom_point(shape = 16, size = 0.3, alpha = 0.5) + 
#  scale_color_gradient(expression(''*log[2]*'(exp)'), low = 'gray', high = 'blue', na.value = 'blue', limits = c(0,3)) + 
#  scale_x_continuous(breaks = seq(-20,20,20)) + scale_y_continuous(breaks = seq(-20,20,20)) + facet_wrap(~variable, ncol = 3) + theme_custom + xlab('tSNE 1') + ylab('tSNE 2')

# Visualization of gene expression: TRC chains
#ts.mk.bd.cxcr3 <- ggplot(df.bd[df.bd$variable %in% 'CXCR3',], aes(x = tx, y = ty, color = value)) + geom_point(shape = 16, size = 0.3, alpha = 0.5) + 
#  scale_color_gradient(expression(''*log[2]*'(exp)'), low = 'gray', high = 'blue', na.value = 'blue', limits = c(0,3)) + 
#  scale_x_continuous(breaks = seq(-20,20,20)) + scale_y_continuous(breaks = seq(-20,20,20)) + facet_wrap(~variable, ncol = 3) + theme_custom + xlab('tSNE 1') + ylab('tSNE 2')

#for(gene in unique(df.bd.all$variable)) {
#ts.all.bd <- ggplot(df.bd.all[df.bd.all$variable == gene, ], aes(x = tx, y = ty, color = value)) + geom_point(shape = 16, size = 0.3, alpha = 0.5) + 
#  scale_color_gradient(low = 'gray', high = 'blue', na.value = 'blue', limits = c(0,3)) + 
#  scale_x_continuous(breaks = c(-25,0,25)) + scale_y_continuous(breaks = c(-25,0,25)) + 
#  facet_wrap(~variable, ncol = 26) + theme_custom + xlab('tSNE 1') + ylab('tSNE 2')
#}
## View plots
#gridExtra::marrangeGrob(grobs = ts.all.bd, nrow = 4, ncol = 5)

# Confidence of annotation
sc <- bd # because function has a problem in calling seurat object, and it calls sc instead of the second argument
SingleR.PlotFeature(sr.bd$singler[[2]]$SingleR.single, bd, plot.feature = 'MaxScore', dot.size = 2)
SingleR.PlotFeature(sr.bd$singler[[2]]$SingleR.single, bd,
                    plot.feature = -log10(sr.bd$singler[[2]]$SingleR.single$pval))
bd.cf <- data.frame(nGene = bd[['nFeature_RNA']],
                    pval = -log10(sr.bd$singler[[2]]$SingleR.single$pval),
                    Identity = sr.bd$singler[[2]]$SingleR.single$labels)
bp.ac.bd <- ggplot(bd.cf, aes(Identity, y = pval, color = Identity)) + geom_boxplot() + scale_color_manual(values = col2) + 
  ylab('-log10(p-value)') + theme(axis.text.x = element_text(angle = 45, hjust = 1)) + theme_custom

# Check number of detected genes per cell
df.bd$ng <- apply(GetAssayData(bd, slot = 'counts'), 2, function(x) sum(x>0))
ts.ng <- ggplot(df.bd, aes(x = tx, y = ty, color = ng)) + 
  geom_point(shape = 16, size = 1, alpha = 0.5) + 
  scale_color_gradientn('', colours = colorRampPalette(rev(RColorBrewer::brewer.pal(n = 11, name = 'RdBu')))(50)[c(1:21,29:50)]) + 
  scale_x_continuous(breaks = c(-25,0,25)) + scale_y_continuous(breaks = c(-25,0,25)) + #facet_wrap(~variable, ncol = 4) + 
  theme_custom + xlab('tSNE 1') + ylab('tSNE 2') + #facet_wrap(~value, ncol = 2) + 
  ggtitle('Number of detected genes per cell')
ts.ng

#-----------------------------------------------------------------------------------------------------------------#
# Differential gene expression visualization (compare to all other cells)
#-----------------------------------------------------------------------------------------------------------------#

# Identify genes differentially expressed between samples (separately) and all other cells
hm <- list()
#genes <- list(sk = list())
DimPlot(bd, reduction = 'tsne', dims = c(1,2), group.by = 'RNA_snn_res.1.2', split.by = 'orig.ident', label = T, label.size = 4, pt.size = 1.5, cols = NULL) + theme_bw() + theme_custom + xlab('tSNE 1') + ylab('tSNE 2') + scale_y_continuous(breaks = seq(-25,25,25)) + scale_x_continuous(breaks = seq(-25,25,25)) #+ NoLegend()
DimPlot(bd, reduction = 'tsne', dims = c(1,2), group.by = 'Sample_Name', split.by = 'orig.ident', label = T, label.size = 4, pt.size = 1.5, cols = NULL) + theme_bw() + theme_custom + xlab('tSNE 1') + ylab('tSNE 2') + scale_y_continuous(breaks = seq(-25,25,25)) + scale_x_continuous(breaks = seq(-25,25,25)) #+ NoLegend()

## Host versus donor cells
### Remove NK and B cells from comparison as they are found only within donor cells
#bd[['group']] <- as.numeric(bd$RNA_snn_res.1.2)
#bd[['group']] <- paste0(stringr::str_split(bd$Sample_Name, "_", n = 2, simplify = T)[,1], '_', bd$RNA_snn_res.1.2)
#bd[['group']] <- paste0(bd$Sample_Name, '_', bd$RNA_snn_res.1.2)
bd[['group']] <- factor(bd$Sample_Name, levels = levels(factor(bd$Sample_Name))[c(1,4,3,6,2,5)] )
#bd[['group']] <- factor(stringr::str_split(bd$Sample_Name, "_", n = 2, simplify = T)[,1], levels = c('blood','skin') )
#bd$group <- plyr::mapvalues(bd$group, from = levels(factor(bd$group)),
#                            to = c('blood', 'blood', 'blood', 'other', 'blood', 'other', 'other', 
#                                   'other', 'other', 'other', 'skin_3', 'other', 'skin_5', 'other'))
#bd$group <- ifelse(bd$group != 'blood' & bd$Sample_Name == 'blood_Th', other, bd$group)
Idents(bd) <- bd$group
#de$skin_3 <- FindMarkers(bd, assay = 'RNA', ident.1 = 'skin_3', ident.2 = 'blood', min.pct = 0.1, 
#                         logfc.threshold = 0.25, test.use = 'wilcox', min.diff.pct = -Inf, random.seed = 1)
#de$skin_3[, 'gene'] <- row.names(de$skin_3)
#de$skin_3 <- (de$skin_3 %>% group_by(avg_logFC<=0) %>% filter(p_val_adj <= 0.05) %>% arrange(p_val) %>% arrange(`avg_logFC <= 0`))
#de$skin_5 <- FindMarkers(bd, assay = 'RNA', ident.1 = 'skin_5', ident.2 = 'blood', min.pct = 0.1, 
#                         logfc.threshold = 0.25, test.use = 'wilcox', min.diff.pct = -Inf, random.seed = 1)
#de$skin_5[, 'gene'] <- row.names(de$skin_5)
#de$skin_5 <- (de$skin_5 %>% group_by(avg_logFC<=0) %>% filter(p_val_adj <= 0.05) %>% arrange(p_val) %>% arrange(`avg_logFC <= 0`))
#de$skin <- FindMarkers(bd, assay = 'RNA', ident.1 = 'skin', ident.2 = 'blood', min.pct = 0.1, 
#                       logfc.threshold = 0.25, test.use = 'wilcox', min.diff.pct = -Inf, random.seed = 1)
#de$skin[, 'gene'] <- row.names(de$skin)
#de$skin <- (de$skin %>% group_by(avg_logFC<=0) %>% filter(p_val_adj <= 0.05) %>% arrange(p_val) %>% arrange(`avg_logFC <= 0`))
st.markers <- FindAllMarkers(bd, min.pct = 0.1, logfc.threshold = 0.25, min.diff.pct = -Inf, only.pos = F, 
                             test.use = 'wilcox', return.thresh = 0.01, max.cells.per.ident = Inf, random.seed = 1)
st.markers <- (st.markers %>% group_by(avg_log2FC<=0) %>% filter(p_val_adj <= 0.05) %>% arrange(p_val) %>% 
                 arrange(cluster) %>% arrange(`avg_log2FC <= 0`))
#Idents(bd) <- bd$RNA_snn_res.1.2

# Export skin markers
#st.markers <- bind_rows(de, .id = 'group')
#write.csv(st.markers, file = paste0(dir.figs, '_tissue_subset_markers.csv'), row.names = T)
st.markers <- read.csv(file = 'GA_AN0289_sc_th_tissue_subset_markers.csv', row.names = 1)

# Visualization of marker genes expression (VlnPlot, FeaturePlot, RidgePlot, CellScatter, and DotPlot)
#genes <- (de$cxcr3 %>% top_n(n = 50, wt = avg_logFC))$gene
#genes$cxcr3 <- (de$cxcr3 %>% group_by(avg_logFC<=0) %>% arrange(-abs(avg_logFC)) %>% arrange(`avg_logFC <= 0`))$gene
#genes$cxcr3 <- (de$cxcr3 %>% group_by(avg_logFC<=0) %>% filter(p_val_adj <= 1) %>% arrange(p_val) %>% arrange(`avg_logFC <= 0`))$gene
#genes$cxcr3.p05 <- (de$cxcr3 %>% group_by(avg_logFC<=0) %>% filter(p_val_adj < 0.05) %>% arrange(p_val) %>% arrange(`avg_logFC <= 0`))$gene
#VlnPlot(bd, features = genes$cxcr3[1:10], cols = NULL, pt.size = 1, adjust = 1,  y.max = NULL, ncol = NULL, slot = 'data')
#FeaturePlot(bd, features = genes$cxcr3[1:10], dims = c(1,2), cols = c('lightgrey', 'blue'), pt.size = NULL, 
#            min.cutoff = NA, max.cutoff = NA, reduction = 'tsne', blend = FALSE, blend.threshold = 0.5, 
#            order = NULL, label = F, label.size = 4, ncol = NULL, coord.fixed = T)
DoHeatmap(bd, features = st.markers[st.markers$avg_logFC>0,]$gene, disp.min = -2.5, disp.max = NULL, slot = 'scale.data', 
          label = F, size = 5.5, hjust = 0.5, angle = 0, combine = T, group.by = 'ident', 
          cells = WhichCells(bd, expression = group %in% c('blood_th1','skin_th1','blood_th2','skin_th2','blood_th17','skin_th17'), invert = F)) + 
  ggtitle('Heatmap of DE genes\n(skin vs. blood Th)')# + NoLegend()

# Organize cells
#cells <- list(skin_3=c(WhichCells(bd, expression = group == "skin_3"), WhichCells(bd, expression = group == "blood")), 
#              skin_5=c(WhichCells(bd, expression = group == "skin_5"), WhichCells(bd, expression = group == "blood")) )
cells <- WhichCells(bd, cells = names(bd$group[order(bd$group)]) )

# Extract scaled data
#df.hm <- list(skin_3 = GetAssayData(bd, slot = "scale.data")[de$skin_3$gene, cells$skin_3], 
#              skin_5 = GetAssayData(bd, slot = "scale.data")[de$skin_5$gene, cells$skin_5])
#df.hm <- GetAssayData(bd, slot = "scale.data")

# Create function to draw annotated heatmap
#col$hm2 <- circlize::colorRamp2(seq(-2, 2, length.out = 40), colorRampPalette(c("magenta", "black", "yellow"))(40))
#col$hm2 <- circlize::colorRamp2(seq(0,4, length.out = 40), viridis::viridis_pal()(40))
col$hm2 <- circlize::colorRamp2(seq(-2, 2, length.out = 40), colorRampPalette(c('royalblue1', 'black', 'yellow'))(40))
draw.heatmap <- function(rq=1, h=2, w=0.5){ 
  #df <- GetAssayData(bd, slot = "scale.data")[st.markers$gene, cells]
  #genes <- c(unique(st.markers[st.markers$avg_logFC>0,]$gene), unique(st.markers[st.markers$avg_logFC<0,]$gene))
  up <- st.markers[st.markers$avg_logFC>0,]
  up <- up[duplicated(up$gene),]
  dn <- st.markers[st.markers$avg_logFC<0,]
  dn <- dn[duplicated(dn$gene),]
  df <- GetAssayData(bd, slot = "scale.data")[rbind(up,dn)$gene, cells]
  #HeatmapAnnotation(genotype = pd$tp.sk$geno.ident[match(cells$hostXdonorCD4.1, names(pd$tp.sk$geno.ident))], col = annColors), 
  #subset_split <- stringr::str_split(subset, "_", n = 2, simplify = T)
  #col.split <- factor(bd$group[match(cells, names(bd$group))], levels = levels(bd$group)[c(1,3,5,2,4,6)])
  col.split <- bd$group[match(cells, names(bd$group))]
  #top.ann <- HeatmapAnnotation(tissue = anno_block(gp = gpar(fill = c("white", "black"), col = 'black', lwd = 1), 
  #                                                 labels = levels(col.split), 
  #                                                 labels_gp = gpar(col = c("black",'white'), fontsize = 12)))
  top.ann <- HeatmapAnnotation(sample = col.split, col = list(sample = setNames(col$ts, levels(bd$group)))) 
  #rowAnnotation(status = annRow$hostXdonorCD4.1[genes$hostXdonorCD4.1,], col = annColors),
  row.split <- factor(ifelse((rbind(up,dn)$avg_logFC > 0) == T, "up", "down"), levels = c("up", "down"))
  left.ann <- rowAnnotation(regulation = anno_block(gp = gpar(fill = c("grey90", "grey90"), lwd = 1),
                                                    labels = unique(row.split), labels_gp = gpar(col = "black", fontsize = 12)))
  title <- paste0('Heatmap of DE genes in Th subsets\nblood and skin')
  hm <- Heatmap(df, cluster_columns = F, cluster_rows = F, name = "z-score",  
                show_row_names = T, show_column_names = F, show_heatmap_legend = T, col = col$hm2, 
                column_split = col.split, row_split = row.split, 
                top_annotation = top.ann, left_annotation = left.ann, column_title = title, 
                column_title_gp = grid::gpar(fontsize = 14), row_title_gp = grid::gpar(fontsize = 0), 
                heatmap_legend_param = list(color_bar = "continuous", title_position = "leftcenter-rot"), 
                row_names_gp = gpar(fontsize = 10), use_raster = T, raster_device = "tiff", raster_quality = rq, 
                height = unit(1*h, "mm")*nrow(df), width = unit(0.1*w, "mm")*ncol(df))
  return(hm)
}
hm$all <- draw.heatmap(1,0.5,0.5)

#-----------------------------------------------------------------------------------------------------------------#
# Differential gene expression visualization (compare to each of the other samples)
#-----------------------------------------------------------------------------------------------------------------#

###############################VENN DIAGRAM here! #####################

#levels(bd$group)[c(1,3,5,2,4,6)]

# Identify genes differentially expressed between samples (separately) and each other sample
de <- list()
for(s1 in levels(bd$group)[c(1,3,5,2,4,6)]) { 
  for(s2 in levels(bd$group)[c(1,3,5,2,4,6)][!(levels(bd$group)[c(1,3,5,2,4,6)] %in% s1)]) { 
   print(paste0('Comparing ', s1, ' vs ', s2))
   de[[paste0(s1,'_vs_',s2)]] <- FindMarkers(bd, assay = 'RNA', ident.1 = s1, ident.2 = s2, min.pct = 0.1, 
        logfc.threshold = 0.25, test.use = 'wilcox', min.diff.pct = -Inf, only.pos = T, random.seed = 1)
   de[[paste0(s1,'_vs_',s2)]][, 'gene'] <- row.names(de[[paste0(s1,'_vs_',s2)]])
   de[[paste0(s1,'_vs_',s2)]][, c('sample', 'ref')]  <- list(s1,s2) }
}
#de[[paste0(s1,'_vs_',s2)]]  <- (de$skin_3 %>% group_by(avg_logFC<=0) %>% filter(p_val_adj <= 0.05) %>% arrange(p_val) %>% arrange(`avg_logFC <= 0`))
des <- bind_rows(de, .id = 'group')
des <- des[des$avg_log2FC > 0 & des$p_val_adj < 0.05,]
des <- des[order(des$gene),]
des <- des[order(des$sample),]
#desf <- plyr::ddply(des, colnames(des)[c(8,7)], summarise, reg='up')
#table(desf$sample_1)

# Find genes which are upregulated in each skin sample compared to all other skin samples
desf <- des[grep('skin', des$sample), ]
#desf <- desf[grep('blood', desf$sample_2), ]
desf <- desf[grep('skin', desf$ref), ]
desf <- plyr::ddply(desf, ~sample+gene, summarize, NA)[,-3]
table(desf$sample)
# Keep only genes with standard genenames
#df.bd3[df.bd3$variable %in% des[des$sample_1 == 'skin_th1',]$gene,]
#desf <- desf[]
write.csv(desf, file = paste0(dir.figs, '_skin_subsets_markers.csv'), row.names = T)
desf <- read.csv(file = paste0(dir.figs, '_skin_subsets_markers.csv'), row.names = 1)
desf

#-----------------------------------------------------------------------------------------------------------------#
# Gene expression in violin plots
#-----------------------------------------------------------------------------------------------------------------#
# Visualization of marker gene expression (violin plot)
df.bd$value.noise <- df.bd$value + rnorm(n = length(df.bd$value)) / 100000 #add noise
#bd.markers$gene[bd.markers$cluster ==3 & bd.markers$p_val_adj < 0.05 & bd.markers$avg_logFC > 0]
vp.bd <- Map(x=unique(df.bd$variable), function(x) {
  qt <- (ceiling(quantile(df.bd[df.bd$variable == x,]$value, 1)) + 
           ceiling(quantile(df.bd[df.bd$variable == x,]$value, 1)) %% 2)
  ggplot(df.bd[df.bd$variable %in% x,], aes(x = cluster, y = value.noise, fill = cluster, color = cluster)) + 
    geom_violin(scale = 'width', size = 0.2, alpha = 0.6, show.legend = F) + 
    geom_jitter(shape = 16, size = 1, alpha = 0.5, position=position_jitter(0.2), color = 'black', show.legend = F) + 
    scale_x_discrete('Clusters', breaks = seq(0,12,1), labels = seq(0,12,1)) + 
    scale_y_continuous(limits = c(-0.1, ifelse(qt>4,qt,4)), breaks = seq(0, ifelse(qt>4,qt,4), length.out = 5)) + 
    scale_fill_manual('Subset', values = col$cluster) + scale_color_manual('Subset', values = col$cluster) + 
    xlab('') + ylab(expression(''*log[2]*'(norm. exp.)')) + 
    ggtitle(paste0('blood and skin Th clusters'), subtitle = 'gene expression') + 
    facet_wrap(~variable, nrow = 3) + theme_custom })
#vp.bd.sig <- Map(x=unique(bd.markers$gene[order(bd.markers$p_val_adj, decreasing = F)]), function(x) {
#  qt <- (ceiling(quantile(df.bd[df.bd$variable == x,]$value, 1)) + 
#           ceiling(quantile(df.bd[df.bd$variable == x,]$value, 1)) %% 2)
#  ggplot(df.bd[df.bd$variable %in% x,], aes(x = cluster, y = value.noise, fill = cluster, color = cluster)) + 
#    geom_violin(scale = 'width', size = 0.2, alpha = 0.6, show.legend = F) + 
#    geom_jitter(shape = 16, size = 1, alpha = 0.5, position=position_jitter(0.2), color = 'black', show.legend = F) + 
#    scale_x_discrete('Clusters', breaks = seq(0,12,1), labels = seq(0,12,1)) + 
#    scale_y_continuous(limits = c(-0.1, ifelse(qt>4,qt,4)), breaks = seq(0, ifelse(qt>4,qt,4), length.out = 5)) + 
#    scale_fill_manual('Subset', values = col$cluster) + scale_color_manual('Subset', values = col$cluster) + 
#    xlab('') + ylab(expression(''*log[2]*'(norm. exp.)')) + 
#    ggtitle(paste0('blood Th clusters'), subtitle = 'gene expression') + 
#    facet_wrap(~variable, nrow = 3) + theme_custom })

vp.bd.ts <- Map(x=unique(df.bd$variable), function(x) {
  qt <- (ceiling(quantile(df.bd[df.bd$variable == x,]$value, 1)) + 
           ceiling(quantile(df.bd[df.bd$variable == x,]$value, 1)) %% 2)
  ggplot(df.bd[df.bd$variable %in% x,], aes(x = tissue_subset, y = value.noise, fill = tissue_subset, color = tissue_subset)) + 
    geom_violin(scale = 'width', size = 0.2, alpha = 0.6, show.legend = T) + 
    geom_jitter(shape = 16, size = 1, alpha = 0.5, position=position_jitter(0.2), color = 'black', show.legend = F) + 
    scale_x_discrete('Sample', breaks = seq(0,12,1), labels = seq(0,12,1)) + 
    scale_y_continuous(limits = c(-0.1, ifelse(qt>4,qt,4)), breaks = seq(0, ifelse(qt>4,qt,4), length.out = 5)) + 
    scale_fill_manual('Sample', values = col$ts) + scale_color_manual('Sample', values = col$ts) + 
    xlab('') + ylab(expression(''*log[2]*'(norm. exp.)')) + 
    ggtitle(paste0('blood and skin Th clusters'), subtitle = 'gene expression') + 
    facet_wrap(~variable, nrow = 3) + theme_custom })

#-----------------------------------------------------------------------------------------------------------------#
# Percentage expressing cells and expression levels per cluster
#-----------------------------------------------------------------------------------------------------------------#

# Prepare data table
df.bd2 <- df.bd[,c('cells','cluster','variable','value')]
df.bd2 <- plyr::ddply(df.bd2, c('cluster','variable'), summarise, 
                      percentage = 100*(sum(value>0)/length(value)), value = mean(value))
df.bd.hi <- plyr::ddply(df.bd2, c('variable'), summarise, max=max(value))
df.bd.hi <- df.bd2[df.bd2$variable %in% df.bd.hi[df.bd.hi$max > 3,]$variable,]
df.bd.lo <- plyr::ddply(df.bd2, c('variable'), summarise, max=max(value))
df.bd.lo <- df.bd2[df.bd2$variable %in% df.bd.lo[df.bd.lo$max <= 3,]$variable,]

# Plot data
#col$hm2 <- colorRampPalette(c('grey80', 'bisque', 'hotpink3', 'royalblue4'))(50)
dot.pct <- list()
dot.pct$bd.hi <- ggplot(df.bd.hi, aes(cluster, variable, color = value, size = percentage)) + 
  geom_point(shape = 16) + 
  scale_x_discrete('clusters', expand = c(0,0.6), position = 'top') + 
  scale_y_discrete('', expand = c(0,0.6), limits = rev(unique(df.bd.hi$variable))) + 
  scale_color_gradientn(expression(''*log[2]*'(exp)'), colours = col$hm, oob = scales::squish,
                        limits = c(2,8), breaks = seq(2,8, length.out = 4)) + 
  scale_size_continuous(limits = c(0,100), breaks = seq(0,100,25), range = c(0,6)) + 
  ggtitle('blood and skin Th clusters', subtitle = 'genes expressed at higher levels') + theme_custom + 
  theme(legend.margin = margin(5,5,5,80), legend.background = element_blank(), 
        plot.margin = unit(c(0.2,0.2,3,0.2),'cm'), aspect.ratio = NULL)
#geom_vline(xintercept = seq(1.5, length(unique(gsea$cont))-0.5, 1), lwd = 0.5, colour = 'grey90') + 
#geom_hline(yintercept = seq(1.5, length(unique(gsea$gs))-0.5, 1), lwd = 0.5, colour = 'grey90') + 
#scale_color_manual('Enriched in:', values = c('firebrick3', 'royalblue'), 
#                   guide = guide_legend(override.aes = list(size = 8))) + 
#scale_size_manual(values = c(4,6,8)) +
#scale_alpha_manual('Adjusted p-value:', values = c(0.1,0.5,1), 
#                   guide = guide_legend(override.aes = list(size = 8))) + 
#scale_fill_manual('Signature:', values = c('white', 'grey70', 'black'), 
#                  guide = guide_legend(override.aes = list(shape = 21, size = 5, 
#                                                           fill = c('white', 'grey70', 'black')))) + 
#coord_flip(xlim = c(1,length(unique(gsea$contrast))), clip = "off") + 
#facet_grid(enriched ~ ., scales = "free", space = "free")
dot.pct$bd.lo <- ggplot(df.bd.lo, aes(cluster, variable, color = value, size = percentage)) + 
  geom_point(shape = 16) + 
  scale_x_discrete('clusters', expand = c(0,0.6), position = 'top') + 
  scale_y_discrete('', expand = c(0,0.6), limits = rev(unique(df.bd.lo$variable))) + 
  scale_color_gradientn(expression(''*log[2]*'(exp)'), colours = col$hm, oob = scales::squish,
                        limits = c(0,3), breaks = seq(0,3, length.out = 4)) + 
  scale_size_continuous(limits = c(0,100), breaks = seq(0,100,25), range = c(0,6)) + 
  ggtitle('blood and skin Th clusters', subtitle = 'genes expressed at lower levels') + theme_custom + 
  theme(legend.margin = margin(5,5,5,80), legend.background = element_blank(), 
        plot.margin = unit(c(0.2,0.2,3,0.2),'cm'), aspect.ratio = NULL)

#-----------------------------------------------------------------------------------------------------------------#
# Percentage expressing cells and expression levels per tissue
#-----------------------------------------------------------------------------------------------------------------#

# Prepare data table
df.bd3 <- data.frame(cells = colnames(bd), 
                     group = bd$group, 
                     t(as.matrix(GetAssayData(bd, slot = 'data'))))#[gn, ])))
df.bd3 <- reshape2::melt(df.bd3, id.vars = c('cells', 'group'))
df.bd3$variable <- sub('\\.','-',df.bd$variable)
df.bd3 <- plyr::ddply(df.bd3, c('group','variable'), summarise, 
                      percentage = 100*(sum(value>0)/length(value)), value = mean(value))
df.bd.hi2 <- plyr::ddply(df.bd3, c('variable'), summarise, max=max(value))
df.bd.hi2 <- df.bd3[df.bd3$variable %in% df.bd.hi2[df.bd.hi2$max > 3,]$variable,]
df.bd.lo2 <- plyr::ddply(df.bd3, c('variable'), summarise, max=max(value))
df.bd.lo2 <- df.bd3[df.bd3$variable %in% df.bd.lo2[df.bd.lo2$max <= 3,]$variable,]

# Plot data
dot.pct$bd.hi2 <- ggplot(df.bd.hi2, aes(group, variable, color = value, size = percentage)) + 
  geom_point(shape = 16) + 
  scale_x_discrete('', expand = c(0,0.6), position = 'top', 
                   limits = unique(df.bd.hi2$group)[!(unique(df.bd.hi2$group) %in% 'other')]) + 
  scale_y_discrete('', expand = c(0,0.6), limits = rev(unique(df.bd.hi2$variable))) + 
  scale_color_gradientn(expression(''*log[2]*'(exp)'), colours = col$hm, oob = scales::squish,
                        limits = c(2,8), breaks = seq(2,8, length.out = 4)) + 
  scale_size_continuous(limits = c(0,100), breaks = seq(0,100,25), range = c(0,6)) + 
  ggtitle('Th clusters', subtitle = 'genes expressed at higher levels') + theme_custom + 
  theme(legend.margin = margin(5,5,5,80), legend.background = element_blank(), 
        plot.margin = unit(c(0.2,0.2,3,0.2),'cm'), aspect.ratio = NULL, 
        axis.text.x = element_text(angle = 45, hjust = 0))
dot.pct$bd.lo2 <- ggplot(df.bd.lo2, aes(group, variable, color = value, size = percentage)) + 
  geom_point(shape = 16) + 
  scale_x_discrete('', expand = c(0,0.6), position = 'top', 
                   limits = unique(df.bd.hi2$group)[!(unique(df.bd.hi2$group) %in% 'other')]) + 
  scale_y_discrete('', expand = c(0,0.6), limits = rev(unique(df.bd.lo2$variable))) + 
  scale_color_gradientn(expression(''*log[2]*'(exp)'), colours = col$hm, oob = scales::squish,
                        limits = c(0,3), breaks = seq(0,3, length.out = 4)) + 
  scale_size_continuous(limits = c(0,100), breaks = seq(0,100,25), range = c(0,6)) + 
  ggtitle('Th clusters', subtitle = 'genes expressed at lower levels') + theme_custom + 
  theme(legend.margin = margin(5,5,5,80), legend.background = element_blank(), 
        plot.margin = unit(c(0.2,0.2,3,0.2),'cm'), aspect.ratio = NULL, 
        axis.text.x = element_text(angle = 45, hjust = 0))

# Plot only significantly dysregulated genes
df.bd.gr <- df.bd3[df.bd3$variable %in% unique(st.markers$gene),]
df.bd.gr.hi <- plyr::ddply(df.bd.gr, c('variable'), summarise, max=max(value))
df.bd.gr.hi <- df.bd.gr[df.bd.gr$variable %in% df.bd.gr.hi[df.bd.gr.hi$max > 3,]$variable,]
df.bd.gr.lo <- plyr::ddply(df.bd.gr, c('variable'), summarise, max=max(value))
df.bd.gr.lo <- df.bd.gr[df.bd.gr$variable %in% df.bd.gr.lo[df.bd.gr.lo$max <= 3,]$variable,]

dot.pct$bd.gr.hi <- ggplot(df.bd.gr.hi, aes(group, variable, color = value, size = percentage)) + 
  geom_point(shape = 16) + 
  scale_x_discrete('', expand = c(0,0.6), position = 'top', 
                   limits = unique(df.bd.hi2$group)[!(unique(df.bd.gr.hi$group) %in% 'other')]) + 
  scale_y_discrete('', expand = c(0,0.6), limits = rev(unique(df.bd.gr.hi$variable))) + 
  scale_color_gradientn(expression(''*log[2]*'(exp)'), colours = col$hm, oob = scales::squish,
                        limits = c(2,8), breaks = seq(2,8, length.out = 4)) + 
  scale_size_continuous(limits = c(0,100), breaks = seq(0,100,25), range = c(0,6)) + 
  ggtitle('Th cells', subtitle = 'DE genes expressed at higher levels') + theme_custom + 
  theme(legend.margin = margin(5,5,5,80), legend.background = element_blank(), 
        plot.margin = unit(c(0.2,0.2,3,0.2),'cm'), aspect.ratio = NULL, 
        axis.text.x = element_text(angle = 45, hjust = 0))
dot.pct$bd.gr.lo <- ggplot(df.bd.gr.lo, aes(group, variable, color = value, size = percentage)) + 
  geom_point(shape = 16) + 
  scale_x_discrete('', expand = c(0,0.6), position = 'top', 
                   limits = unique(df.bd.hi2$group)[!(unique(df.bd.gr.lo$group) %in% 'other')]) + 
  scale_y_discrete('', expand = c(0,0.6), limits = rev(unique(df.bd.gr.lo$variable))) + 
  scale_color_gradientn(expression(''*log[2]*'(exp)'), colours = col$hm, oob = scales::squish,
                        limits = c(0,3), breaks = seq(0,3, length.out = 4)) + 
  scale_size_continuous(limits = c(0,100), breaks = seq(0,100,25), range = c(0,6)) + 
  ggtitle('Th cells', subtitle = 'DE genes expressed at lower levels') + theme_custom + 
  theme(legend.margin = margin(5,5,5,80), legend.background = element_blank(), 
        plot.margin = unit(c(0.2,0.2,3,0.2),'cm'), aspect.ratio = NULL, 
        axis.text.x = element_text(angle = 45, hjust = 0))

#-----------------------------------------------------------------------------------------------------------------#
# Save results
#-----------------------------------------------------------------------------------------------------------------#

# Function to save plots with same size
grid <- function(data, layout=rbind(c(1,1,1,1,2,3,3,3))) {
  grid.arrange(grobs = list(data + guides(color = 'none', fill = 'none', size = 'none'), 
                            cowplot::get_legend(data)), ncol = 2, layout_matrix = layout, 
               top = grid::textGrob('', gp = grid::gpar(fontsize = 16, font = 2)))
}
sv.list <- list(ts.cl = grid(ts.cl), um.cl = grid(um.cl), 
                ts.ct = grid(ts.ct), um.ct = grid(um.ct), 
                dot.bd.hi = grid(dot.pct$bd.hi, rbind(c(1,1,1,1,2,3))), 
                dot.bd.lo = grid(dot.pct$bd.lo, rbind(c(1,1,1,1,2,3))), 
                dot.bd.hi2 = grid(dot.pct$bd.hi2, rbind(c(1,1,1,1,2,3))), 
                dot.bd.lo2 = grid(dot.pct$bd.lo2, rbind(c(1,1,1,1,2,3))), 
                dot.bd.gr.hi = grid(dot.pct$bd.gr.hi, rbind(c(1,1,1,1,2,3,3))), 
                dot.bd.gr.lo = grid(dot.pct$bd.gr.lo, rbind(c(1,1,1,1,2,3))))
#dot.sk.hi = grid(dot.pct$sk.hi, rbind(c(1,1,1,1,2,3,3))), 
#dot.sk.lo = grid(dot.pct$sk.lo, rbind(c(1,1,1,1,2,3,3))))

# Save plots
ggsave(plot = sv.list$ts.cl, file = paste0(dir.figs, '_tsne_clusters.pdf'), height = 4, width = 8)
ggsave(plot = sv.list$um.cl, file = paste0(dir.figs, '_umap_clusters.pdf'), height = 4, width = 8)
ggsave(plot = sv.list$ts.ct, file = paste0(dir.figs, '_tsne_tissue.pdf'), height = 4, width = 8)
ggsave(plot = sv.list$um.ct, file = paste0(dir.figs, '_umap_tissue.pdf'), height = 4, width = 8)
ggsave(plot = sv.list$dot.bd.hi, file = paste0(dir.figs, '_pct_expression_clusters_high.pdf'), height = 12, width = 5.6)
ggsave(plot = sv.list$dot.bd.lo, file = paste0(dir.figs, '_pct_expression_clusters_low.pdf'), height = 68, 
       width = 6.6, limitsize = F)
ggsave(plot = sv.list$dot.bd.hi2, file = paste0(dir.figs, '_pct_expression_tissue_high.pdf'), height = 8.6, width = 3.6)
ggsave(plot = sv.list$dot.bd.lo2, file = paste0(dir.figs, '_pct_expression_tissue_low.pdf'), height = 73, 
       width = 4.7, limitsize = F)
ggsave(plot = sv.list$dot.bd.gr.hi, file = paste0(dir.figs, '_pct_expression_de_high.pdf'), height = 7.6, width = 4.2)
ggsave(plot = sv.list$dot.bd.gr.lo, file = paste0(dir.figs, '_pct_expression_de_low.pdf'), height = 23, width = 4.5)

ggsave(paste0(dir.figs, '_tsne_expression.pdf'), 
       gridExtra::marrangeGrob(grobs = ts.all.bd, nrow = 3, ncol = 2), height = 11, width = 10)
ggsave(paste0(dir.figs, '_umap_expression.pdf'), 
       gridExtra::marrangeGrob(grobs = um.all.bd, nrow = 3, ncol = 2), height = 11, width = 10)
ggsave(paste0(dir.figs, '_violin_expression.pdf'), 
       gridExtra::marrangeGrob(grobs = vp.bd, nrow = 3, ncol = 2), height = 11, width = 10)
ggsave(paste0(dir.figs, '_violin_expression_samples.pdf'), 
       gridExtra::marrangeGrob(grobs = vp.bd.ts, nrow = 3, ncol = 2), height = 11, width = 10)

# Save complex heatmaps
pdf(paste0(dir.figs, "_heatmap_samples.pdf"), height = 15, width = 12)
draw(draw.heatmap(rq=1, h=3, w=0.5))
dev.off()

# Save plots for QC analysis
ggsave(plot = vp.bdc, file = paste0(dir.figs.qc, '_violin_expression_cluster.pdf'), height = 5, width = 15)
ggsave(plot = vf.bd, file = paste0(dir.figs.qc, '_var_features.pdf'), height = 4, width = 6)
ggsave(plot = hm.pc.bd, file = paste0(dir.figs.qc, '_heatmap_pca.pdf'), height = 36, width = 20)
ggsave(plot = hm.pc.bd2, file = paste0(dir.figs.qc, '_heatmap_pca2.pdf'), height = 28, width = 20)
ggsave(plot = ep.bd, file = paste0(dir.figs.qc, '_elbowplot.pdf'), height = 4, width = 4)
ggsave(plot = hm.mk.bd, file = paste0(dir.figs.qc, '_heatmap_markers.pdf'), height = 12, width = 20)
#ggsave(plot = bp.ac.bd, file = paste0(dir.figs.qc, '_boxplot_ann_confidence.pdf'), height = 5, width = 5)

#tt <- ttheme_minimal(rowhead=list(fg_params=list(fontface='bold')))
#pdf(paste0(dir.figs, '_cell_number.pdf'), height = 18, width = 10)
#grid.arrange(tableGrob(sr.bdn1, theme=tt), tableGrob(sr.bdn2, theme=tt), 
#             tableGrob(sr.bdn3, theme=tt), tableGrob(sr.bdn4, theme=tt), nrow = 4, ncol = 1)
#dev.off()

#-----------------------------------------------------------------------------------------------------------------#
# End of script
#-----------------------------------------------------------------------------------------------------------------#
