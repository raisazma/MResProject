library(dplyr)
library(Seurat)
library(patchwork)
library(R.utils)
library(ggplot2)

cerebellum.data <- readRDS("~/cerebellardata/GSM5952337_Fetal_cerebellum_F (1).rds")

DimPlot(cerebellum.data, reduction = "umap", label = TRUE, group.by = "cluster.names")

# QC
cerebellum.data[["percent.mt"]] <- PercentageFeatureSet(cerebellum.data, pattern = "^MT-")

# Show QC metrics for the first 5 cells
head(cerebellum.data@meta.data, 5)

#View(cerebellum.data@meta.data)

# Visualize QC metrics as a violin plot
VlnPlot(cerebellum.data, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3)

# FeatureScatter 
FeatureScatter(cerebellum.data, feature1 = "nCount_RNA", feature2 = "nFeature_RNA") + geom_smooth(method = 'lm')
plot1 <- FeatureScatter(cerebellum.data, feature1 = "nCount_RNA", feature2 = "percent.mt")
plot2 <- FeatureScatter(cerebellum.data, feature1 = "nCount_RNA", feature2 = "nFeature_RNA")
plot1 + plot2

cerebellum.data <- subset(cerebellum.data, subset = nFeature_RNA > 200 &  percent.mt < 10)

# Normalize data
cerebellum.data <- NormalizeData(cerebellum.data)

# plot variable features with and without labels
cerebellum.data <- FindVariableFeatures(cerebellum.data, selection.method = "vst", nfeatures = 2000)

# Identify the 10 most highly variable genes
top10 <- head(VariableFeatures(cerebellum.data), 10)

# plot variable features with and without labels
plot1 <- VariableFeaturePlot(cerebellum.data)
plot2 <- LabelPoints(plot = plot1, points = top10, repel = TRUE)
plot1 + plot2

all.genes <- rownames(cerebellum.data)
cerebellum.data <- ScaleData(cerebellum.data, features = all.genes)

# perform linear dimensional reduction
cerebellum.data <- RunPCA(cerebellum.data, features = VariableFeatures(object = cerebellum.data))

# Examine and visualize PCA results a few different ways
print(cerebellum.data[["pca"]], dims = 1:5, nfeatures = 5)

VizDimLoadings(cerebellum.data, dims = 1:2, reduction = "pca")
DimPlot(cerebellum.data, reduction = "pca")
DimHeatmap(cerebellum.data, dims = 1:9, cells = 500, balanced = TRUE)

ElbowPlot(cerebellum.data)
cerebellum.data <- FindNeighbors(cerebellum.data, dims = 1:15)
cerebellum.data <- FindClusters(cerebellum.data, resolution = 0.5)

# Look at cluster IDs of the first 5 cells
head(Idents(cerebellum.data), 5)

#### store the current clustering in a variable Idents(cerebellum.data) -> cerebellum.data$my_new_clustering

# run UMAP
cerebellum.data <- RunUMAP(cerebellum.data, dims = 1:10)
DimPlot(cerebellum.data, reduction = "umap", label = TRUE)
DimPlot(cerebellum.data, reduction = "umap", label = TRUE, group.by = "cluster.names")

# test out feature and violin plot
FeaturePlot(cerebellum.data, "MYC")
VlnPlot(cerebellum.data, features = c("SOX2","SOX11"), group.by = "cluster.names")

### give me the names of all the clusters
levels(cerebellum.data$cluster.names)

### these are the categories I want to retain
clusters.of.interest <- c("01NSC", "02TCP", "06GCP_pro", "07GCP1","08GCP2","09Granule","10UBC_PRO","11UBC","12UBCdiff")

### subset the columns of this large seurat object to give only the cluster we want
subset.cerebellum.data <- cerebellum.data[,cerebellum.data$cluster.names%in%clusters.of.interest]

### these are the genes which define 02 RL precursors by aldinger
VlnPlot(subset.cerebellum.data, features = c("OTX2","LMX1A","MKI67"), group.by = "cluster.names")

VlnPlot(subset.cerebellum.data, features = c("OTX2","LMX1A","GRIA1","PTPRK","MEIS2","RBFOX3","ITPR1","DNAH6","SOX2","DNAH11","ITPK1"), group.by = "cluster.names")

### TCP markers as defined by Lu Lab paper 02TCP
VlnPlot(subset.cerebellum.data,features = c("HNRNPH1","SOX11"), group.by = "cluster.names")

### genes listed from Supp Fig2 panel D in Richard Lu's, including some used to define the aldinger set
VlnPlot(subset.cerebellum.data, features = c("WLS","LMX1A","OTX2","HNRNPH1","SOX11","RSPO3","EOMES","BARHL1"),group.by = "cluster.names")

### genes listed from Supp Fig2 panel D in Richard Lu's, including some used to define the aldinger set
VlnPlot(subset.cerebellum.data, features = c("SOX2","NMU","CIP2A","KPNA2","UNCX"),group.by = "cluster.names")
#FeaturePlot(subset.cerebellum.data, "OTX2")                                            

### extract the values for the projection
Assays(subset.cerebellum.data, "RNA")

### extract the normalised count data
tpms.mat <- subset.cerebellum.data[["RNA"]]@scale.data

# install and load libraries

install.packages('NMF', dependencies = TRUE)
install.packages("MASS")
install.packages("NMF")

# load required data objects
nmf.res <- readRDS(file = "/home/raisaa/Group3-4App-main/StarProtocols_Guide/data/nmf.res.rds")
source(file = "/home/raisaa/Group3-4App-main/StarProtocols_Guide/R/Project_NMF.R")

# load sample data as matrix object

# project NMF model to sequencing data
tpms.mat <- project.NMF(input.array = as.matrix(tpms.mat), nmf.result = nmf.res)

# rows 3 and 1 in tpms.mat correspond to the Groups 3 and 4 metagenes, extract them from the data and transpose matrix
g3g4.tpms <- t(tpms.mat[c(3,1),]) 

# apply logistic transformation to metagenes
logistic.g3g4.tpms <- apply(g3g4.tpms,2,function(x){(1 / (1 + exp(-x)))}) 
logistic.g3g4.tpms.score <- apply(logistic.g3g4.tpms,1,function(x){x[2]/(x[1]+x[2])})

# scale values between 0 to 1
logistic.g3g4.tpms.continuum.score <- (logistic.g3g4.tpms.score-min(logistic.g3g4.tpms.score)) / (max(logistic.g3g4.tpms.score)-min(logistic.g3g4.tpms.score))
logistic.g3g4.tpms.continuum.score <- logistic.g3g4.tpms.score

View(logistic.g3g4.tpms.continuum.score)

# present output as data frame for export
logistic.g3g4.tpms.continuum.score <- as.data.frame(logistic.g3g4.tpms.continuum.score)
colnames(logistic.g3g4.tpms.continuum.score)
data.frame(logistic.g3g4.tpms.continuum.score)

# projected results, map them back onto the seurat object
logistic.g3g4.tpms.continuum.score <- subset.cerebellum.data$logistic.g3g4.tpms.continuum.score
View(logistic.g3g4.tpms.continuum.score) ## each row should represent sample, column corresponds to continuum score value

as.numeric(tpms.mat[1,]) -> subset.cerebellum.data$group3_metagene
as.numeric(tpms.mat[3,]) -> subset.cerebellum.data$group4_metagene
as.numeric(tpms.mat[4,]) -> subset.cerebellum.data$shh_metagene
as.numeric(tpms.mat[2,]) -> subset.cerebellum.data$wnt_metagene
logistic.g3g4.tpms.continuum.score -> subset.cerebellum.data$g3g4.score

VlnPlot(subset.cerebellum.data, features = "group4_metagene", group.by = "cluster.names")
VlnPlot(subset.cerebellum.data, features = "group3_metagene" ,group.by = "cluster.names")
VlnPlot(subset.cerebellum.data, features = "shh_metagene" ,group.by = "cluster.names")
VlnPlot(subset.cerebellum.data, features = "wnt_metagene" ,group.by = "cluster.names")
VlnPlot(subset.cerebellum.data, features = "wnt_metagene")
VlnPlot(subset.cerebellum.data, features = "shh_metagene")
VlnPlot(subset.cerebellum.data, features = "group3_metagene")
VlnPlot(subset.cerebellum.data, features = "group4_metagene")

### draw some dimplots that show the clusters and the metagene and the g3g4 score ###confused
DimPlot(object = subset.cerebellum.data, group.by = "cluster.names") -> dimplot.cluster.names
FeaturePlot(object = subset.cerebellum.data,  "group4_metagene", cols = c("grey","green"), min.cutoff = 0) -> dimplot.group4
FeaturePlot(object = subset.cerebellum.data, "group3_metagene", cols = c("grey","yellow"), min.cutoff = 0) -> dimplot.group3
FeaturePlot(object = subset.cerebellum.data, "shh_metagene", cols = c("grey","red"),min.cutoff=0) -> dimplot.shh
FeaturePlot(object = subset.cerebellum.data, "wnt_metagene", cols = c("grey","blue"), min.cutoff=0) -> dimplot.wnt

dimplot.cluster.names
dimplot.group4
dimplot.group3
dimplot.shh
dimplot.wnt

library(gridExtra)
grid.arrange(dimplot.cluster.names,
             dimplot.group4,
             dimplot.group3,
             dimplot.shh,
             dimplot.wnt,
             nrow = 3
) -> dimplot.all

dimplot.all

ggsave(dimplot.cluster.names, file = "./dimplot.cluster.names.png")
ggsave(dimplot.group4, file = "./dimplot.group4.png")
ggsave(dimplot.group3, file = "./dimplot.group3.png")
ggsave(dimplot.shh, file = "./dimplot.wnt.png")
ggsave(dimplot.wnt, file = "./dimplot.shh.png")
ggsave(dimplot.all, file = "./dimplot.all.png", height = 12, width = 8)

DimPlot(object = subset.cerebellum.data$wnt_metagene, group.by = "cluster.names")

VizDimLoadings(subset.cerebellum.data, dims = 1:2, reduction = "pca")
DimPlot(subset.cerebellum.data$group4_metagene, reduction = "pca", group.by = "cluster.names")
DimPlot(subset.cerebellum.data$group3_metagene, group.by = "cluster.names")

# Examine PCA for the new subset
subset.cerebellum.data <- RunPCA(subset.cerebellum.data, features = VariableFeatures(object = subset.cerebellum.data))
print(subset.cerebellum.data[["pca"]], dims = 1:5, nfeatures = 5)

### Visualize PCA results
VizDimLoadings(subset.cerebellum.data, dims = 1:2, reduction = "pca")
DimPlot(subset.cerebellum.data, reduction = "pca", group.by = "cluster.names")
DimHeatmap(subset.cerebellum.data, dims = 1:9, cells = 500, balanced = TRUE)
ElbowPlot(subset.cerebellum.data)
subset.cerebellum.data <- JackStraw(subset.cerebellum.data, num.replicate = 100)
subset.cerebellum.data <- ScoreJackStraw(subset.cerebellum.data, dims = 1:20)

#### Eecalculate the UMAP with the subset and also recalculate the clusters and the PCA and the variable genes...

subset.cerebellum.data <- RunUMAP(subset.cerebellum.data, dims = 1:10)
DimPlot(subset.cerebellum.data, reduction = "umap", label = TRUE, group.by = "cluster.names")
subset.cerebellum.data <- FindVariableFeatures(subset.cerebellum.data, selection.method = "vst", nfeatures = 2000)

### make a new clustering and directly compare your clustering with the Lu paper clustering
subset.cerebellum.data <- FindNeighbors(subset.cerebellum.data, dims = 1:15)
subset.cerebellum.data <- FindClusters(subset.cerebellum.data, resolution = 0.5)
View(subset.cerebellum.data)

##### Ident(subset.cerebellum.data) will be the new clustering that you made
View(subset.cerebellum.data$cluster.names)
#Idents(subset.cerebellum.data) <- subset.cerebellum.data$cluster.names
head(Idents(subset.cerebellum.data), 5)
View(Idents(subset.cerebellum.data))

### subset.cerebellum.data$cluster.names is the name of there clustering
### you can use the function table(a,b) to create a comparison table
### prop.table(a,b)

table(Idents(subset.cerebellum.data), subset.cerebellum.data$cluster.names)
prop.table(table(Idents(subset.cerebellum.data), (subset.cerebellum.data$cluster.names)), 2) -> proport

library("pheatmap")

pheatmap(proport, cluster_rows = F, cluster_cols = F, scale = "none")

