#load package into R
library(immunarch)

#path to files
file_path <- c("./Data/10XData_Skin/", "./Data/10XData_Blood/")

names(file_path)
# Output: NULL

names(file_path) <- sapply(1:length(file_path), function(i) paste0("Sample", i))
names(file_path)
# Output: Sample1 Sample2 ... Sample10

#load 10x
immdata <- repLoad(file_path)

#look at data
immdata

#exploring dataset
exp_len <- repExplore(immdata$data, .method = "len", .col = "aa")
exp_cnt <- repExplore(immdata$data, .method = "count")
exp_vol <- repExplore(immdata$data, .method = "volume")

p1 <- vis(exp_len)
p2 <- vis(exp_cnt)
p3 <- vis(exp_vol)

p1
p2 + p3

#get most abundant clonotypes
top(immdata$data[[1]])

#filter clonotypes
#coding(immdata$data[[1]])
#noncoding(immdata$data[[1]])
#nrow(inframes(immdata$data[[1]]))
#nrow(outofframes(immdata$data[[1]]))

#downsampling
#ds = repSample(immdata$data, "downsample", 100)
#sapply(ds, nrow)
#ds = repSample(immdata$data, "sample", .n = 10)
#sapply(ds, nrow)


# 1. Analyse
exp_len <- repExplore(immdata$data, .method = "len", .col = "aa")
# 2. Visualise
p1 <- vis(exp_len)
# 3. Fix and make publication-ready results
fixVis(p1)

#Clonality
imm_pr <- repClonality(immdata$data, .method = "clonal.prop")
imm_pr

#highly prolific
imm_top <- repClonality(immdata$data, .method = "top", .head = c(10, 100, 1000, 3000, 10000))
imm_top

#least prolific
imm_rare <- repClonality(immdata$data, .method = "rare")
imm_rare

vis(imm_top) + vis(imm_rare)

#homeostasis
imm_hom <- repClonality(immdata$data,
                        .method = "homeo",
                        .clone.types = c(Small = .0001, Medium = .001, Large = .01, Hyperexpanded = 1)
)
imm_hom
vis(imm_hom)

#repertoire overlap
imm_ov1 <- repOverlap(immdata$data, .method = "public", .verbose = F)
imm_ov2 <- repOverlap(immdata$data, .method = "morisita", .verbose = F)
imm_ov3 <- repOverlap(immdata$data, .method = "jaccard", .verbose = F)

p1 <- vis(imm_ov1)
p2 <- vis(imm_ov2, .text.size = 2)
p3 <- vis(imm_ov3)

p1
p2
p3
#p1 + p2

vis(imm_ov1, "heatmap2")
p1 <- vis(imm_ov2, .text.size = 2.5, .signif.digits = 1)
p2 <- vis(imm_ov2, .text.size = 2, .signif.digits = 2)

p1 + p2


# Apply different analysis algorithms to the matrix of public clonotypes:
# "mds" - Multi-dimensional Scaling
repOverlapAnalysis(imm_ov1, "mds")



# Pass "nt" as the second parameter to build the public repertoire table using CDR3 nucleotide sequences
pr.nt <- pubRep(immdata$data, "nt", .verbose = F)
pr.nt


# Pass "aa+v" as the second parameter to build the public repertoire table using CDR3 aminoacid sequences and V alleles
# In order to use only CDR3 aminoacid sequences, just pass "aa"
pr.aav <- pubRep(immdata$data, "aa+v", .verbose = F)
pr.aav


#Repertoire diversity
# Compute statistics and visualise them
# gini diversity measure
div_gini <- repDiversity(immdata$data, "gini.simp")

# Hill numbers
div_hill <- repDiversity(immdata$data, "hill")

# D50
div_d50 <- repDiversity(immdata$data, "d50")

# Ecological diversity measure
div_div <- repDiversity(immdata$data, "div")

p1 <- vis(div_gini)
p3 <- vis(div_hill)
p4 <- vis(div_d50)
p6 <- vis(div_div)

p1
p3
p4
p6

# Load the dplyr library
library(dplyr)

# Load the database with immunarch
vdjdb = dbLoad("/Users/mahim/Desktop/HiWi/RScripts/Immunoarch/Data/Databases/VDJDB/SearchTable-2021-07-08 18_28_16.96.tsv", "vdjdb-search", .species = "HomoSapiens", .chain = "TRB")
vdjdb

mcpas = dbLoad("/Users/mahim/Desktop/HiWi/RScripts/Immunoarch/Data/Databases/McPAS-TCR/McPAS-TCR.csv", "mcpas", .species = "Human", .chain = "TRB", .pathology = "Epstein Barr virus (EBV)")
mcpas

tbadb = dbLoad("/Users/mahim/Desktop/HiWi/RScripts/Immunoarch/Data/Databases/TBAdb/TBAdb.xlsx", "tbadb", .species = "Homo Sapiens", .chain = c("TRB", "TRA-TRB"), .pathology = "EBV")
tbadb

# Check which antigen epitopes are presented in the database
#epitope <- vdjdb$Epitope
#vis(epitope, .title = "Epitopes")
table(vdjdb$Epitope)
table(mcpas$Epitope.peptide)



#Repertoire annotaion of found clones
dbAnnotate(immdata$data, vdjdb, "CDR3.aa", "CDR3")
dbAnnotate(immdata$data, mcpas, c("CDR3.aa", "V.name"), c("CDR3.beta.aa", "TRBV"))
imm_gu <- dbAnnotate(immdata$data, vdjdb, "CDR3.aa", "CDR3")
#vis_hist(imm_gu, .title = "Clonotype Annotation")

