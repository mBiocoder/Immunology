#=================================================================================================================#
# Exploratory analysis of BD Rhapsody targeted scRNA-Seq from Th subsets (CD25-, CD45RA-, CCR6, CXCR3, CCR4) 
# isolated from healthy skin and blood
# Date: 11.11.2021
# Author: Mahima Arunkumar (most of the code adapted from Gustavo)
# Adapted from https://satijalab.org/seurat/pbmc3k_tutorial.html
# Adapted from https://www.biorxiv.org/content/biorxiv/suppl/2018/03/19/284604.DC1/284604-2.pdf
# Rscript purpose:
# - analyze BD Rhapsody data processed with Seven Bridges Genomics platform for the two replicates
# - gene set enrichment analysis using Seurat scoring for Trm and circulating signatures
#Human CD4 CD103 cutaneous resident memory T cells are found in the circulation of healthy individuals, Maria M. Klicznik 2019
#=================================================================================================================#
#=================================================================================================================#
# Prepare R session
#=================================================================================================================#

# Load libraries
#.libPaths(paste0(.libPaths(), '/Seurat_v3.1.1'))#, paste0(.libPaths(), '/Seurat_v3.1.1'))
library(openxlsx)
library(Seurat)
library(ggplot2)
library(dplyr)
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
experimentid2 <- c('EX0014')

# Define Analysis_ID according to ngs_analysis_list
analysisid <- sub('.R', '', stringr::str_split(rstudioapi::getSourceEditorContext()$path, '/', simplify = T)[,5])

# Define samples analyzed (all)
#sa <- c('cd')

# Define reference sample name
#ct <- c('blood')
#tr <- c('skin')

# Define |log2fc| cutoff
lfc.cutoff <- 0

# Define analysis description
an.desc <- c('sc_bd_thsubsets_trm_gratz')
an.descs <- c('sc_th_trm_gratz')

#-----------------------------------------------------------------------------------------------------------------#
# Load samples and analysis data
#-----------------------------------------------------------------------------------------------------------------#

# Load table containing the ngs sample list filtered for the data analyzed here
nsl <- read.xlsx('./ngs_sample_list.xlsx', sheet = 1, startRow = 1, colNames = T,
                 rowNames = F, skipEmptyRows = T, skipEmptyCols = T, rows = NULL, cols = NULL, 
                 check.names = T, namedRegion = NULL, na.strings = 'NA', fillMergedCells = F)
nsl <- nsl[nsl$Experiment_ID == experimentid1 | nsl$Experiment_ID == experimentid2, ]
nsl <- plyr::ddply(nsl, colnames(nsl)[1:18], plyr::summarise, sum = NA)

# Update table containing the analyses list to include this analysis
# Include info like which read alignment algorithm was used
# Load, update avoiding duplicates (when script is run several times) and save back overwriting old one
nal <- read.xlsx('./ngs_analysis_list.xlsx', sheet = 1, startRow = 1, colNames = T,
                 rowNames = F, skipEmptyRows = T, skipEmptyCols = T, rows = NULL, cols = NULL, 
                 check.names = T, namedRegion = NULL, na.strings = 'NA', fillMergedCells = F)
nal <- nal[nal$Analysis_ID == analysisid & !is.na(nal$Analysis_ID), ]

# Sample data
metadata <- data.frame(row.names = unique(nsl$Sample_name), 
                       reshape2::colsplit(unique(nsl$Sample_name), '_', c('exp', 'type', 'subset')))
#head(metadata)

#-----------------------------------------------------------------------------------------------------------------#
# Archive old runs of this same analysis
#-----------------------------------------------------------------------------------------------------------------#


# Create new folder to store figures from this analysis
dir.create(paste0('figures/', analysisid, '_', an.desc))
dir.create(paste0('figures/', analysisid, '_', an.desc, '/rep1'))
dir.create(paste0('figures/', analysisid, '_', an.desc, '/rep2'))
dir.figs <- list(all = paste0('figures/', analysisid, '_', an.desc, '/', analysisid, '_', an.descs) )
dir.figs$rep1 <- paste0('figures/', analysisid, '_', an.desc, '/rep1/', analysisid, '_', an.descs, '_rep1')
dir.figs$rep2 <- paste0('figures/', analysisid, '_', an.desc, '/rep2/', analysisid, '_', an.descs, '_rep2')
#dir.figs.pca <- paste0('figures/', analysisid, '_', an.desc, '/pca/', analysisid, '_', an.descs)

#=================================================================================================================#
# Analysis of processed data
#=================================================================================================================#
#-----------------------------------------------------------------------------------------------------------------#
# Subset distribution
#-----------------------------------------------------------------------------------------------------------------#

# Load Seurat objects
bd <- list(#rep1=readRDS(paste0('../201911/figures/GA_AN0289_sc_bd_blood_skin_th/GA_AN0289_sc_th_seurat.rds')),
  rep1=readRDS(paste0('./GA_AN0316_sc_th_rep1_seurat.rds')), 
  rep2=readRDS(paste0('./GA_AN0310_sc_th_rep2_seurat.rds')))
#bd$rep1$subset <- factor(bd$rep1$Sample_Name[drop = T], levels = levels(bd$rep1$Sample_Name[drop = T])[c(1,3,2,4,6,5)])
#bd$rep2$subset <- factor(bd$rep2$Sample_Name[drop = T], levels = levels(bd$rep2$Sample_Name[drop = T])[c(1,3,2,4,6,5)])
bd$rep1$subset <- plyr::mapvalues(bd$rep1$Sample_Name[drop = T], from = levels(factor(bd$rep1$Sample_Name))[c(1,3,2,4,6,5)],
                                  to = c('bTh1','bTh17','bTh2','sTh1','sTh17','sTh2'))
#bd$rep1$subset <- factor(bd$rep1$subset, levels = levels(bd$rep1$subset)[c(1,3,2,4,6,5)])
#bd$rep1$subset <- factor(bd$rep1$subset, levels = levels(bd$rep1$subset)[c(1,4,3,6,2,5)])
bd$rep2$subset <- plyr::mapvalues(bd$rep2$Sample_Name[drop = T], from = levels(factor(bd$rep2$Sample_Name))[c(1,3,2,4,6,5)],
                                  to = c('bTh1','bTh17','bTh2','sTh1','sTh17','sTh2'))
#bd$rep2$subset <- factor(bd$rep2$subset, levels = levels(bd$rep2$subset)[c(1,3,2,4,6,5)])
#bd$rep2$subset <- factor(bd$rep2$subset, levels = levels(bd$rep2$subset)[c(1,4,3,6,2,5)])

# Change cluster numbers
bd$rep1$cluster <- plyr::mapvalues(bd$rep1$seurat_clusters, from = levels(bd$rep1$seurat_clusters), to = c(1:8))
bd$rep2$cluster <- factor(plyr::mapvalues(bd$rep2$seurat_clusters[drop=T], from = levels(bd$rep1$seurat_clusters), to = c(2,1,4,6,3,5,7,8)), levels=1:8)

# Create palette of colors
col <- list(subset = c(RColorBrewer::brewer.pal(n = 12, name = 'Paired')[c(1:12)], 'grey60'), 
            hm = colorRampPalette(c('grey80', 'skyblue2', 'royalblue1', 'royalblue4'))(50), 
            hm2 = colorRampPalette(c('royalblue1', 'black', 'yellow'))(40))
#col
#=================================================================================================================#
# GSEA
#=================================================================================================================#
#-----------------------------------------------------------------------------------------------------------------#
# Prepare gene list
#-----------------------------------------------------------------------------------------------------------------#

# Gene list mackay
trm <- c('SEMA6D','ITGAE','IFITM2','CDH1','GZMB','ST6GAL1','NR4A3','HSPA1A','LGALS3','LMNA','GEM','GCNT1','ANXA1','CD244','AREG','STYK1','GSTM3','XCL2','DENND4A','RGS2','AHR','RBPJ','RGS1','CAPG','KRT5','HAVCR2','PTPN5','QPCT','LITAF','TJP1','ANXA2','C15orf48','GZMA','EGR2','BCL2A1','INPP4B','NR4A1','ANXA5','ADGRG1','TNF','CXCR6','VIM','ITGA1','VPS37B','DENND5A','CLNK','CD44','FNDC3A',
         'ATP5IF1','CD3G','RGCC','TIPARP','GPR171','ADAM8','GADD45B','FCER1G','CISH','ICOS','TNFSF9','JUP','NFIB','HSP90AA1','REL','PDCD1','CDKN1A','SWAP70','EGR1','FGL2','CSRNP1','PYGL','DAPK2','NR4A2','MRNIP','PTPRJ','MYADM','MED10','ARNT2','DKKL1','SRGN','ADGRG3','FOSL2','SIK1B','CEP85L','IFNG','OSBPL3','C21orf91','TOB2','GSTK1','SH3BGRL3','ANAPC1','RASGEF1B','CTNND2','NUDT4B','NEURL3',
         'CMKLR1','HSPA5','ABCB1','SYTL3','TOP2A','CTSW','S100A4','LIN54','BAIAP2','JUNB','CD6','PHLDA1','TMEM154','KDM6B','IRF4','CSF1','CTLA4','BTG2','PRKCH','GPR141','ENDOD1','TRIB1','DDX3X','IGFLR1','ISG20','FBXO30','HASPIN','GCH1','ARHGAP21','SLC16A10','ZFAND5','PPP1CC','SERPINB9','TUBB2A','IFRD1','PLK3','ORAI1','TNFRSF1B','CTNNA1','SQSTM1','KRT15','IL2RA','LY6G5B','GPR34','NFKBID','LDLRAD4','DNAJB6','HSP90AB1','CD96','DUSP1','SLC39A6','ELL2','UBC','LPAR6','TRAV23DV6','ADORA3','DNAJA1','GBF1','ODC1','SRSF3','CCL4','NDFIP1','CCL18','SLC38A2','TMPO','GALNT3','SOAT1','TP53INP1','CTSB','S100A11','IL4R','PRKRIP1','SKIL','MT2A','RHOB','RAB7A','PTP4A1','NFE2L2','EMP1','ZNF318','FKBP2','ISY1-RAB43','HMGN5','PHACTR2','B4GALNT4','PRMT2','HSPA8','CLCN3','HEXIM1','ARL5B','ITSN1','WSB1','RDM1','AC134669.1','TPPP3','SLK','PTPRE','SLC6A6','ERF','SLAMF7','UBL3','FOS','TNFRSF9','EID1','PELO','RIN3','TMEM30A','RELL1','PERP','DLD','KCTD12','ZFP36','SCD','PTK2B','KLF4','TGFBR3','CSTB','ECHDC2','TOPORS','IL2RB','TES','WDR78','NABP1','LAT2','S100A6','TNFSF10','PRNP','VBP1','DSG1','UBE2J1','CD14','JMY','PPP2R2D','PANX1','SNX9','TMEM123','SLC41A2','P2RY10','AGAP1','COPG1','TIAM1','IER2','PTGFRN','MYO6','PPP1R3B','NFIL3','CNTN1','FAM174B','PTMA','NT5E','IER5','KIAA1109','PDIA6','TRPM6','STK17B','SERTAD1','HSPH1','HSPB1','TIGIT','PRKAR1A','PBRM1','DNAJB4','PLD3','SARAF','CD226','WNK1','DNAJB1','COQ10B','KIF13B','P2RY14','BRD2','PCNA','SLBP','CLIC1','EIF4E3','RNF149','UBB','MYO1E','CALR','EYA2','HSPA13','ZBTB20','KPNA3','FASLG','SPTY2D1')

#original mackay
cir <- c('SELL','S1PR1','ITGA4','GIMAP4','CD55','SLAMF6','GIMAP7','SIDT1','CCL5','ADGRE5','EOMES','ITGB1','PLAC8','USP33','GMFG','GIMAP1-GIMAP5','LEF1','DGKA','RASGRP2','KLF3','ACP5','TLR1','S1PR5','GRAMD4','GIMAP6','ATIC','IL7R','DNAH8','SLFN12','ZNF845','PRKCQ','STK38','PYHIN1','TCF7','LPIN1','CNOT6L','AC008878.3','TMEM131L','MS4A6A','TXK','XRN2','ARHGAP4','PODNL1','SLFN5','NSG2',
         'NUP210','EPSTI1','KLHL6','MDN1','NOD1','RASA3','CD84','PDE2A','SMPDL3B','KLRC4-KLRK1','INIP','PPP1R12A','RNF138','SAMHD1','GAB3','RPL13','PDK1','NFATC3','UTRN','MACF1','GPR18','RPS5','CD22','SP110','KBTBD11','DENND2D','CD27','GZMM','ZC4H2','C1orf21','KLRG1','ELOVL7','GSTP1','RPS16','DPP4','GSAP','TBC1D1','TRAF3IP3','LYST','DERL1','UBAP2L','BICDL1','RPL32','CCDC91','MGST2',
         'HERC6','CPSF7','RPL18A','RAB5B','LIMD2','RPS2','DAPL1','RIPOR2','KIF23','AC104389.5','SRCAP','XDH','CD2AP','BACH2','B3GNT5','ABHD11','SUPT5H','SATB1','APBB1IP','RACK1','RPL38','PLCB2','GIMAP8','IL27RA','ASAP1','EVL','MRPL2','SH3BP5','EGFL8','BRWD1','ZNF397','CD79A','ZNF429','ATP10D','STAT4','ZC3H12D','RPL29','WSCD2','THTPA','AFF3','FAM160A2','C9orf85','PEX11G','TAF1D','RPL18',
         'IL18RAP','PPDPF','SPTBN1','SAMD3')

#genes discussed with Christina
#genes <- c('S1PR1', 'SELL', 'CD69', 'CD103', 'LGALS3', 'RUNX3', 'FABP4', 'FABP5')

#Bd_rhapsody gene panel
#data <- c('FCGR1A','FCGR1B','FCGR1C','CD68','LAMP1','P2RX7','PANX1','SLC2A1','SLC2A2','SLC2A3','SLC2A4','SLC5A12','CD36','SLC27A1','SLC27A2','SLC27A4','FABP1','FABP2','FABP3','FABP4','FABP5','FABP6','FABP7','SLC29A1','SLC29A3','SLC7A5','SLC7A1','SLC1A5','CD14','FAS','CD19','TCF7','ABCC11')

#958_zielisnki genepanel
#cir <- c('ABCB1','ABCC1','AHR','ANTXR1','AOC3','AREG','ATXN1','BACH2','BATF','BCL2L1','BIRC3','CAPG','CASP1','CASP12','CASP8','CCL17','CCL22','CD101','CD28','CD34','CD38','CRTAM','CYP1A1','DEPTOR','FABP5','FUT7','GBP1','GCLC','HAVCR1','HIF1A','HSPA1A','ID2','ID3','IFNAR1','IFNAR2','IGF1','IKBKB','IKBKE','IKZF3','IL10','IL1A','IL1B','IRAK4','ITGA1','ITGB1','ITGB7','JAK1','JAML','KLF2','KRT14','MAF','MAP3K7','MAPK1','MET','MS4A1','MTOR','MYD88','NCAM1','NFAT5','NLRP3','NOD1','NOD2','NOTCH1','NR4A1','NRP1','P2RX7','PAG1','PRKCD','PTK7','PYCARD','RARA','RELA','RIPK1','RPS4Y1','RPTOR','RXRA','S100A4','SGK1','SLC2A1','STAM2','STAT2','TAB1','TBK1','TLR4','TLR7','TLR8','TMEM173','TOX','TRAF2','TRAF5','TSLP','TYK2','TYR','VDAC1','XCL2','XIST','ZBTB7B','ZEB2')

#human t cell panel
#cir <- c('AIM2','ANXA5','APOBEC3G','ARL4C','AURKB','B3GAT1','BAX','BCL11B','BCL2','BCL6','BIN2','BTG1','BTLA','C10orf54','CASP3','CBLB','CCL1','CCL2','CCL20','CCL3','CCL4','CCL5','CCNB1','CCND2','CCR1','CCR10','CCR2','CCR3','CCR4','CCR5','CCR6','CCR7','CCR8','CCR9','CD160','CD2','CD244','CD247','CD27','CD274','CD300A','CD3D','CD3E','CD3G','CD4','CD40LG','CD44','CD48','CD5','CD52','CD6','CD69','CD7','CD70','CD8A','CD8B','CD9','CHI3L2','CLC','CLEC2D','CNOT2','CSF2','CSF3','CST7','CTLA4','CTSW','CX3CR1','CXCL10','CXCL13','CXCL8','CXCL9','CXCR1','CXCR3','CXCR4','CXCR5','CXCR6','DPP4','DUSP1','DUSP2','DUSP4','EGR1','EGR3','ENTPD1','EOMES','F5','FAS','FASLG','FBXO22','FOSB','FOSL1','FOXO1','FOXO3','FOXP1','FOXP3','FYB','FYN','GAPDH','GATA3','GHR','GIMAP2','GIMAP5','GIMAP7','GLG1','GNLY','GZMA','GZMB','GZMH','GZMK','GZMM','HAVCR2','HLA-A','HLA-C','HLA-DMA','HLA-DMB','HLA-DPA1','HLA-DPB1','HLA-DQA1','HLA-DQB1','HLA-DRA','HLA-DRB3','HMGB2','HMMR','ICAM1','ICOS','IER3','IER5','IFNG','IFNGR1','IKZF2','IL12A','IL12RB1','IL12RB2','IL13','IL15','IL15RA','IL17A','IL17F','IL18','IL18R1','IL18RAP','IL1R2','IL2','IL21','IL22','IL23R','IL25','IL2RA','IL2RB','IL3','IL31','IL32','IL4','IL4R','IL5','IL6','IL6R','IL7R','IL9','IL9R','IRF4','IRF8','ITGA4','ITGAE','ITGAL','ITGAM','ITGAX','ITGB2','ITK','JUN','JUNB','KIT','KLRB1','KLRC1','KLRC3','KLRC4','KLRF1','KLRG1','KLRK1','LAG3','LAIR2','LAP3','LAT','LAT2','LCK','LEF1','LGALS1','LGALS3','LIF','LILRB4','LRRC32','LTA','LTB','MKI67','MYC','NAMPT','NCR3','NINJ2','NKG7','NT5E','OAS1','PASK','PDCD1','PECAM1','PIK3IP1','PMCH','PRDM1','PRF1','PTGDR2','PTTG2','PYCR1','RORA','RORC','RUNX3','S1PR1','SELL','SELPLG','SEMA7A','SLAMF1','SPOCK2','SPP1','STAT1','STAT3','STAT4','STAT5A','STAT6','TARP_refseq','TBX21','TCF7','TGFB1','TGFB3','TIAF1','TIGIT','TK1','TLR2','TLR9','TNF','TNFRSF18','TNFRSF1B','TNFRSF25','TNFRSF4','TNFRSF8','TNFRSF9','TNFSF10','TOP2A','TRAC','TRAT1','TRBC2','TRDC','TRIB2','TSPAN32','TXK','TYMS','UBE2C','VNN2','XCL1','ZAP70','ZBED2','ZBTB16','ZNF683')


#trm[!(trm %in% rownames(bd$rep1))] <- c('TENT5A')
#trm <- trm[(trm %in% rownames(bd$rep1))]
#cir[!(cir %in% rownames(bd$rep1))] <- c('TENT5A')
#cir <- cir[(cir %in% rownames(bd$rep1))]
#gs <- list(trm = trm[(trm %in% rownames(bd$rep1))], cir = cir[(cir %in% rownames(bd$rep1))])

# Define gene sets
set.seed(1)
gs <- list(trm = sort(as.character(trm)),
           cir = sort(as.character(cir)), 
           random = sample(rownames(bd$rep1), 15))
rm(.Random.seed)
gs$trm <- gs$trm[gs$trm %in% rownames(bd$rep1)]
gs$cir <- gs$cir[gs$cir %in% rownames(bd$rep1)]

gs$trm
gs$cir


#-----------------------------------------------------------------------------------------------------------------#
# Calculate module scores using Seurat
#-----------------------------------------------------------------------------------------------------------------#

# Calculate module scores (remove genes not detected in any cell for binning to work)
bd <- Map(x=names(bd), function(x) { 
  all <- AddModuleScore(object = subset(bd[[x]], features = rownames(bd)[
    rowSums(as.matrix(GetAssayData(bd[[x]], slot = 'data'))) > 0]), features = gs, 
    ctrl = ceiling(0.2*median(apply(GetAssayData(bd[[x]], slot = 'counts'), 2, function(z) sum(z>0)))), nbin = 24, seed = 1, 
    name = 'scores')
  AddMetaData(object = bd[[x]], metadata = all@meta.data[,grep('scores',colnames(all@meta.data))], col.name = names(gs)) })


# Create data table
#for(i in names(gs)) { 
#df.ms[[i]] <- Map(x=names(bd), function(x) {
#  df <- reshape2::melt(data.frame(tx = bd[[x]]@reductions$tsne@cell.embeddings[,1], 
#                                  ty = bd[[x]]@reductions$tsne@cell.embeddings[,2], 
#                                  ux = bd[[x]]@reductions$umap@cell.embeddings[,1], 
#                                  uy = bd[[x]]@reductions$umap@cell.embeddings[,2], 
#                                  sample = bd[[x]]@meta.data$sample, 
#                                  cells = colnames(bd[[x]]), 
#                                  cluster = bd[[x]]@meta.data$cluster,
#                                  subset = bd[[x]]@meta.data$subset, 
#                                  celltype = bd[[x]]@meta.data$celltype, 
#                                  genotype = bd[[x]]@meta.data$genotype, 
#                                  #setNames(as.data.frame(GetAssayData(bd$bd, slot = "data")[c('RPS4Y1'),]), nm = 'RPS4Y1')) )
#                                  t(as.matrix(GetAssayData(bd[[x]], slot = 'data')[gs[[i]],]))), 
#                       id.vars = c('tx', 'ty', 'ux', 'uy', 'sample', 'cells', 'cluster', 'subset', 'celltype', 'genotype')) 
#  cbind(df, value.noise = (df$value + rnorm(n = length(df$value)) / 100000)) })
#}

# Plot scores
ts.ms <- list()
#ts.ms
for(i in names(gs)) { 
  ts.ms[[i]] <- Map(x=names(bd), function(x) { 
    df <- data.frame(tx = bd[[x]]@reductions$tsne@cell.embeddings[,1], ty = bd[[x]]@reductions$tsne@cell.embeddings[,2], 
                     #trm_scores = scales::rescale(bd[[x]]@meta.data$trm_score1,to=c(0,1)), 
                     #circulating_scores = scales::rescale(bd[[x]]@meta.data$circulating_score2,to=c(0,1)))
                     score = scales::rescale(bd[[x]]@meta.data[[i]],to=c(0,1)))
    ggplot(df, aes(x = tx, y = ty, color = score)) + 
      geom_point(shape = 16, size = 1, alpha = 1) + 
      #scale_color_gradient(low = 'gray', high = 'blue', na.value = 'blue', limits = c(min(df.mk[[k]][x,]$value),3)) + 
      scale_color_gradientn('Score', colours = col$hm, oob = scales::squish, 
                            limits = c(0,1), breaks = seq(0,1, length.out = 3)) + 
      scale_x_continuous(breaks = seq(-20,20,20)) + scale_y_continuous(breaks = seq(-25,25,25)) + 
      #facet_wrap(~variable+donor, ncol = 4) + 
      theme_custom + xlab('tSNE 1') + ylab('tSNE 2') + 
      ggtitle(paste0('blood and skin Th subsets'), subtitle = paste0(i)) })
}

#-----------------------------------------------------------------------------------------------------------------#
# Compare gene set scores from different subsets and do t-test 
#-----------------------------------------------------------------------------------------------------------------#
#Uncomment the following chunks of code 
# Compare enrichment within subsets to get the respecttive p-values

#Comparison of bTh1-bTh17
#ttest <- Map(x=names(bd), function(x) { lapply(setNames(names(gs), nm = names(gs)), function(y) { 
  #sapply(setNames(c('bTh1-bTh17'), nm = c('bTh1-bTh17')), function(z) { 
    #formatC(t.test(bd[[x]][[y]][bd[[x]]$subset=="bTh1",], bd[[x]][[y]][bd[[x]]$subset=="bTh17",])$p.value, 
            #format = 'g', digits = 2)
  #})})})

#print all p-values for the boxplot 
#ttest

#Comparison of bTh17-bTh2
#ttest <- Map(x=names(bd), function(x) { lapply(setNames(names(gs), nm = names(gs)), function(y) { 
  #sapply(setNames(c('bTh17-bTh2'), nm = c('bTh17-bTh2')), function(z) { 
    #formatC(t.test(bd[[x]][[y]][bd[[x]]$subset=="bTh17",], bd[[x]][[y]][bd[[x]]$subset=="bTh2",])$p.value, 
            #format = 'g', digits = 2)
  #})})})

#print all p-values for the boxplot 
#ttest

#Comparison of bTh1-bTh2
#ttest <- Map(x=names(bd), function(x) { lapply(setNames(names(gs), nm = names(gs)), function(y) { 
#sapply(setNames(c('bTh1-bTh2'), nm = c('bTh1-bTh2')), function(z) { 
#formatC(t.test(bd[[x]][[y]][bd[[x]]$subset=="bTh1",], bd[[x]][[y]][bd[[x]]$subset=="bTh2",])$p.value, 
#format = 'g', digits = 2)
#})})})

#print all p-values for the boxplot 
#ttest

#Comparison of sTh1-sTh17
#ttest <- Map(x=names(bd), function(x) { lapply(setNames(names(gs), nm = names(gs)), function(y) { 
#sapply(setNames(c('sTh1-sTh17'), nm = c('sTh1-sTh17')), function(z) { 
#formatC(t.test(bd[[x]][[y]][bd[[x]]$subset=="sTh1",], bd[[x]][[y]][bd[[x]]$subset=="sTh17",])$p.value, 
#format = 'g', digits = 2)
#})})})

#print all p-values for the boxplot 
#ttest

#Comparison of sTh17-sTh2
#ttest <- Map(x=names(bd), function(x) { lapply(setNames(names(gs), nm = names(gs)), function(y) { 
#sapply(setNames(c('sTh17-bTh2'), nm = c('sTh17-sTh2')), function(z) { 
#formatC(t.test(bd[[x]][[y]][bd[[x]]$subset=="sTh17",], bd[[x]][[y]][bd[[x]]$subset=="sTh2",])$p.value, 
#format = 'g', digits = 2)
#})})})

#print all p-values for the boxplot 
#ttest

#Comparison of sTh1-sTh2
ttest <- Map(x=names(bd), function(x) { lapply(setNames(names(gs), nm = names(gs)), function(y) { 
sapply(setNames(c('sTh1-sTh2'), nm = c('sTh1-sTh2')), function(z) { 
formatC(t.test(bd[[x]][[y]][bd[[x]]$subset=="sTh1",], bd[[x]][[y]][bd[[x]]$subset=="sTh2",])$p.value, 
format = 'g', digits = 2)
})})})

#print all p-values for the boxplot 
ttest


# Box plot for scores from donor and host cells
df.ss <- histplot.ss <- boxplot.ss <- list()
for(geneSetName in names(gs)) {
  df.ss[[geneSetName]] <-  Map(x=names(bd), function(x) { 
    data.frame(tx = bd[[x]]@reductions$tsne@cell.embeddings[,1], ty = bd[[x]]@reductions$tsne@cell.embeddings[,2], 
               cell = names(Idents(bd[[x]])), #genotype = bd[[x]]$genotype, 
               subset = bd[[x]]$subset, 
               #subset <- factor(bd[[x]]$subset, levels = c("bTh1", "bTh2", "bTh17", "sTh1", "sTh2", "sTh17")),
               score = scales::rescale(bd[[x]]@meta.data[[geneSetName]],to=c(0,1))) })
  boxplot.ss[[geneSetName]] <- Map(x=names(bd), function(x) { 
    ggplot(df.ss[[geneSetName]][[x]], aes_string(x='subset', y='score', color = 'subset', fill = 'subset')) + 
      geom_jitter(size = 1, shape = 21, color = 'grey40', alpha = 0.6, position=position_jitter(0.2), show.legend = F) + 
      geom_boxplot(outlier.colour = NA, color = 'grey40', width=0.6, alpha = 0.4, show.legend = F) + 
      #geom_text(data=, aes(y=0.95, label=ttest[[x]][[geneSetName]])) + 
      #geom_text(data = dmf2, aes(y = 0.95), #position = position_dodge(width = 0.7), 
      #          nudge_x = 0, size = 8, color = "black") + #position = position_dodge(width = 0.7), 
      annotate("text", label = paste0('p = ',ttest[[x]][[geneSetName]]), x = c(1.5,3.5,5.5), y = 0.88, size = 3.8, colour = "black") + 
      #geom_segment(data = dmf2, colour = "black", show.legend = F, 
      #             aes(x = ifelse(sig == "*", as.numeric(genename) - 0.2, NA), y = (mean+se + 0.3), 
      #                 xend = ifelse(sig == "*", as.numeric(genename) + 0.2, NA), yend = (mean+se + 0.3))) + 
      annotate("segment", x = c(1,3,5), xend = c(2,4,6), y = 0.85, yend = 0.85, size = 0.5, colour = "black") + 
      scale_color_manual('', values = col$subset) + #scale_shape_manual('', values = c(21,24)) + 
      scale_fill_manual('', values = col$subset) + 
      #guides(shape = guide_legend(override.aes = list(size = 4, alpha = 1))) + 
      scale_x_discrete('') + scale_y_continuous('Score', breaks = seq(0,1,0.2), limits = c(-0.01,1.01)) + 
      ggtitle(paste0('blood and skin Th subsets'), subtitle = paste0(gsub('_',' ',geneSetName))) + theme_custom + 
      theme(aspect.ratio = NULL, axis.text.x = element_text(size = 12.8, color = "black", angle = 45, hjust = 1))#+ facet_wrap(~subset, ncol = 5, scales = 'fixed')
  })
  
  # Plot frequency of scores for host and donor cells
  histplot.ss[[geneSetName]] <- Map(x=names(bd), function(x) { 
    max <- ceiling(max(sapply(levels(df.ss[[geneSetName]][[x]]$subset), function(y) max(density(df.ss[[geneSetName]][[x]][df.ss[[geneSetName]][[x]]$subset %in% y,'score'])$y))))
    ggplot(df.ss[[geneSetName]][[x]], aes_string(x='score', color='subset', fill='subset')) + 
      #geom_histogram(color = 'black', breaks = seq(0,0.25,0.05), position = 'identity', alpha = 0.3, size = 0.2) + 
      geom_density(aes(y=..density..), position = 'identity', alpha = 0.3, size = 0.2) + 
      #scale_x_continuous('Number of cells containing a SNP', sec.axis = dup_axis(name = '', breaks = median(rowSums(alt.sk>0)), labels = round(median(rowSums(alt.sk>0))))) + #, limits = c(-0.05,1.05), breaks = seq(0,1,0.2), labels = function(n){format(n, scientific = F)}, expand = c(0,0)) + 
      scale_x_continuous('Score', limits = c(0,1), breaks = c(seq(0,1,0.5))) + 
      #scale_y_continuous('Density', limits = c(0,6), breaks = c(seq(0,6,6))) + 
      scale_y_continuous('Density', limits = c(0,max), breaks = c(0,max)) + 
      scale_fill_manual('Subset', values = col$subset) + 
      scale_color_manual('Subset', values = col$subset) + 
      ggtitle(paste0('blood and skin Th subsets'), subtitle = paste0(gsub('_',' ',geneSetName))) + theme_custom + 
      theme(aspect.ratio = NULL, strip.text = element_text(size = 16, face = 'bold', margin = margin(0,0,0,5))) + 
      facet_wrap(~subset, strip.position = 'right', ncol = 1, scales = 'fixed')
  })
}

boxplot.ss

#Heatmap

# Store plots
#hist[[i]] <- recordPlot()
#heat[[i]] <- pheatmap::pheatmap(cells.assign[[i]], cluster_rows = T, show_rownames = T, cluster_cols = T, show_colnames = F, 
                                #border_color = NA, #annotation_col = ann_col, annotation_colors = ann_colors, 
                                #col = c('grey90', 'black'), legend = F, 
                                #legend_breaks = seq(0,1,1), main = 'assigned cells', silent = T)
#heat[[i]]

head(trm)

DoHeatmap(bd, features = trm, disp.min = -2.5, disp.max = NULL, slot = 'scale.data', 
          label = F, size = 5.5, hjust = 0.5, angle = 0, combine = T, group.by = 'ident', 
          cells = WhichCells(bd, expression = group %in% c('blood_th1','skin_th1','blood_th2','skin_th2','blood_th17','skin_th17'), invert = F)) + 
  ggtitle('Heatmap of DE genes\n(skin vs. blood Th)')# + NoLegend()



#gn <- c('IL7R', 'CCR7', 'CD4', 'CD8A', 'CD8B', 'GNLY', 'NKG7', 'GZMB', 'TRDC')
# Visualization of marker genes expression (VlnPlot, FeaturePlot, RidgePlot, CellScatter, and DotPlot)
#VlnPlot(bd, features = bd.top$gene[1:6], cols = NULL, pt.size = 1, adjust = 1,  y.max = NULL, ncol = NULL, slot = 'data')
#VlnPlot(bd, features = gn, cols = NULL, pt.size = 1, adjust = 1,  y.max = NULL, ncol = NULL, slot = 'data')

#siehe auch andere zeilen in GA_AN0316 ab zeile 290
#Sonst try seaborn.violinplot und maybe auch heatmap, da es für object of class "list" auch gehen sollte


#-----------------------------------------------------------------------------------------------------------------#
# Save results
#-----------------------------------------------------------------------------------------------------------------#

# Save  tSNE plots
lapply(names(ts.ms), function(x) lapply(names(ts.ms[[x]]), function(y) 
  ggsave(paste0(dir.figs[[y]], '_tsne_score_', x, '.pdf'), ts.ms[[x]][[y]], height = 4, width = 4.5)) )

# Save boxplots
lapply(names(boxplot.ss), function(x) lapply(names(boxplot.ss[[x]]), function(y) 
  ggsave(paste0(dir.figs[[y]], '_boxplot_subset_', x, '.pdf'), boxplot.ss[[x]][[y]], height = 4, width = 4) ))

# Save histograms
lapply(names(boxplot.ss), function(x) lapply(names(boxplot.ss[[x]]), function(y) 
  ggsave(paste0(dir.figs[[y]], '_hist_subset_', x, '.pdf'), histplot.ss[[x]][[y]], height = 6, width = 4.5) ))

#-----------------------------------------------------------------------------------------------------------------#
# End of script