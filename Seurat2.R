# loading Seurat
library(Seurat)
library(tidyverse)

#importing data
NML1 <- Read10X(data.dir = "GSE132771_RAW/NML1/")
NML2 <- Read10X(data.dir = "GSE132771_RAW/NML2/")
NML3 <- Read10X(data.dir = "GSE132771_RAW/NML3/")

# Create of Seurat object
NML1 <- CreateSeuratObject(counts = NML1, project = "NML1", min.cells = 3, min.features = 200)
NML2 <- CreateSeuratObject(counts = NML2, project = "NML2", min.cells = 3, min.features = 200)
NML3 <- CreateSeuratObject(counts = NML3, project = "NML3", min.cells = 3, min.features = 200)

# View our Seurat objects
NML1
rownames(NML1)
colnames(NML1)
view(NML1)
view(NML1@meta.data)

# save Seurat object
saveRDS(NML1, file ="/Users/rezamoosavi_1/Desktop/RNAseq/Seurat Practice 2/NML1.RDS")
saveRDS(NML2, file ="/Users/rezamoosavi_1/Desktop/RNAseq/Seurat Practice 2/NML2.RDS")
saveRDS(NML3, file ="/Users/rezamoosavi_1/Desktop/RNAseq/Seurat Practice 2/NML3.RDS")

# read RDS file again into R
NML1 <- readRDS("NML1.RDS")
NML2 <- readRDS("NML2.RDS")
NML3 <- readRDS("NML3.RDS")

# Merge Seurat objects
mergedNML <- merge(NML1, y = c(NML2, NML3),
                   add.cell.ids = ls()[1:3],
                   project = "MergedNML")
ls()
mergedNML
view(mergedNML@meta.data)

# Save merged Seurat object
saveRDS(mergedNML, file = "/Users/rezamoosavi_1/Desktop/RNAseq/Seurat Practice 2/mergedNML.RDS")

# Quality control features (n_features, n_counts, percentage_mt)
view(mergedNML@meta.data)
range(mergedNML$nFeature_RNA)
range(mergedNML$nCount_RNA)

# store mitocondrial percentage in object metadata
mergedNML <- PercentageFeatureSet(mergedNML, pattern = "^MT-", col.name = "percent.mt")
view(mergedNML@meta.data)
range(mergedNML$percent.mt)

# now we can use a violin plot to check the distribution of our columns in the metadata
VlnPlot(mergedNML, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3)

# check the number of cells currently exist in our object
mergedNML


# now let's get rid of low quality cells
mergedNML <- subset(mergedNML, subset = nFeature_RNA > 500 & nFeature_RNA <4000 & nCount_RNA < 20000 & percent.mt < 10)

VlnPlot(mergedNML, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3)

mergedNML

# now we need to normalize our data.
view(mergedNML)
mergedNML@assays[["RNA"]]@data@x

mergedNML <- NormalizeData(mergedNML, normalization.method = "LogNormalize", scale.factor = 10000)

mergedNML@assays[["RNA"]]@data@x

# Calculate the set of highly variable features
mergedNML
mergedNML <- FindVariableFeatures(mergedNML, selection.method = "vst", nfeatures = 2000)

# features and features plots
mergedNML@assays[["RNA"]]@var.features
VariableFeaturePlot(mergedNML)

# Data Scaling
mergedNML <- ScaleData(mergedNML) # based on 2000 found variable genes

all.genes <- rownames(mergedNML)
mergedNML <- ScaleData(mergedNML, features = all.genes)

mergedNML@assays[["RNA"]]@scale.data
mergedNML@commands

# save the Seurat object up to this point
saveRDS(mergedNML, file = "/Users/rezamoosavi_1/Desktop/RNAseq/Seurat Practice 2/mergedNML.RDS")

#performing PCA
mergedNML <- readRDS("mergedNML.RDS")
mergedNML

mergedNML <- RunPCA(mergedNML)

# now lets plot our data with dimplot
# there are also vizDimReduction() and DimHeatmap () to do this in Seurat
DimPlot(mergedNML, reduction = "pca", dims = c(1, 2))
DimPlot(mergedNML, reduction = "pca", dims = c(1, 10))
DimPlot(mergedNML, reduction = "pca", dims = c(1, 50))

# determining PCAs with greatest explained variance (like scree plot)

# the second method is called ElbowPlot function
ElbowPlot(mergedNML)
ElbowPlot(mergedNML, ndims = 50, reduction = "pca")

# there are several methods to do this
#mergedNML <- JackStraw(mergedNML, num.replicate = 100)
#mergedNML <- ScoreJackStraw(mergedNML, dims = 1:20)
#JackStrawPlot(mergedNML, dims = 1:20)

# Clustering cells
mergedNML <- FindNeighbors(mergedNML, dims = 1:20)
mergedNML <- FindClusters(mergedNML, resolution = 0.1)
mergedNML <- FindClusters(mergedNML, resolution = 0.3)

# Using non-linear dimentionality reduction using tSNE/UMAP
mergedNML <- RunUMAP(mergedNML, dims = 1:20)
DimPlot(mergedNML, reduction = "umap", label = TRUE, repel = TRUE)

# now let's cluster cells with tSNE
mergedNML <- RunTSNE(object = mergedNML)
DimPlot(object = mergedNML, reduction = "tsne")

# now let's check the structure of our seurat object
view(mergedNML)

# save for further analysis
saveRDS(mergedNML, file = "/Users/rezamoosavi_1/Desktop/RNAseq/Seurat Practice 2/mergedNML.RDS")
mergedNML <- readRDS("mergedNML.RDS")




