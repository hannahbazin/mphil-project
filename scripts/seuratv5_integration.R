# This code attempts to do integration on three samples using the IntegrateLayers() function from Seurat v5.
# It does not separate the samples appropriately, as they mostly cluster by sample of origin.
# Therefore I use Seurat v4 for integration (see file human_PVAT_snRNAseq_analysis.Rmd).

# Load libraries
library(GEOquery)
library(Seurat)
library(dplyr)
library(patchwork)
library(ggplot2)
library(ggpubr)
library(tidyverse)

# Set working directory: setwd("/Users/hannahbazin/Desktop/Cambridge/Academics/Han_Lab/MPhil/mphil-project")

# Load data
human1_data <- Read10X(data.dir = "data/GSE164528_RAW/GSM5068996_hPVAT1/")
human2_data <- Read10X(data.dir = "data/GSE164528_RAW/GSM5068997_hPVAT2/")
human3_data <- Read10X(data.dir = "data/GSE164528_RAW/GSM5068998_hPVAT3/")

# Create Seurat objects
human1 <- CreateSeuratObject(counts = human1_data, project = "GSE164528", min.cells = 3, min.features = 200)
human2 <- CreateSeuratObject(counts = human2_data, project = "GSE164528", min.cells = 3, min.features = 200)
human3 <- CreateSeuratObject(counts = human3_data, project = "GSE164528", min.cells = 3, min.features = 200)

# See number of cells in each sample
ncol(human1)
ncol(human2)
ncol(human3)

# Pre-processing workflow
human1[["percent.mt"]] <- PercentageFeatureSet(human1, pattern = "^MT-")
human2[["percent.mt"]] <- PercentageFeatureSet(human2, pattern = "^MT-")
human3[["percent.mt"]] <- PercentageFeatureSet(human3, pattern = "^MT-")

VlnPlot(human1, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3)
VlnPlot(human2, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3)
VlnPlot(human3, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3)

plot1 <- FeatureScatter(human1, feature1 = "nCount_RNA", feature2 = "percent.mt")
plot2 <- FeatureScatter(human1, feature1 = "nCount_RNA", feature2 = "nFeature_RNA")
plot1 + plot2

plot1 <- FeatureScatter(human2, feature1 = "nCount_RNA", feature2 = "percent.mt")
plot2 <- FeatureScatter(human2, feature1 = "nCount_RNA", feature2 = "nFeature_RNA")
plot1 + plot2

plot1 <- FeatureScatter(human3, feature1 = "nCount_RNA", feature2 = "percent.mt")
plot2 <- FeatureScatter(human3, feature1 = "nCount_RNA", feature2 = "nFeature_RNA")
plot1 + plot2

human1 <- subset(human1, subset = nFeature_RNA > 1000 & nFeature_RNA < 6000 & nCount_RNA > 500 & nCount_RNA < 40000 & percent.mt < 1)
human2 <- subset(human2, subset = nFeature_RNA > 1000 & nFeature_RNA < 6000 & nCount_RNA > 500 & nCount_RNA < 40000 & percent.mt < 1)
human3 <- subset(human3, subset = nFeature_RNA > 1000 & nFeature_RNA < 6000 & nCount_RNA > 500 & nCount_RNA < 40000 & percent.mt < 1)

# Seurat v5 integration
# Step 0: Normalise and identify variable features for each dataset
human1 <- human1 %>% NormalizeData() %>% 
  FindVariableFeatures()
human2 <- human2 %>% NormalizeData() %>% 
  FindVariableFeatures()
human3 <- human3 %>% NormalizeData() %>% 
  FindVariableFeatures()

# Step 1: Create a list of Seurat objects
ser_list <- list(human1, human2, human3)

# Step 2: Merge the Seurat objects into one
ser_merged <- merge(
  ser_list[[1]],
  y = c(ser_list[[2]], ser_list[[3]]),
  add.cell.ids = c("human1", "human2", "human3"),
  project = "GSE164528"
)

# Step 3: Perform standard normalization and pre-integration analysis
ser_pre <- ser_merged %>%
  NormalizeData() %>%
  FindVariableFeatures() %>%
  # add this to ScaleData() later (takes a long time to run):
  # vars.to.regress = 'percent.mt'
  ScaleData() %>%
  RunPCA(., features = VariableFeatures(object = ser_merged)) %>%
  FindNeighbors(dims = 1:10) %>%
  FindClusters(resolution = 0.5) %>%
  RunUMAP(dims = 1:10, reduction = "pca", reduction.name = "umap.unintegrated")


# Step 4: Perform integration using layers (new method in Seurat v5.1.0)
ser_int <- IntegrateLayers(
  object = ser_pre,
  method = CCAIntegration,
  orig.reduction = "pca",
  new.reduction = "integrated.cca",
  verbose = FALSE
)

# Step 5: Finalize integration
ser_int <- ser_int %>%
  FindNeighbors(reduction = "integrated.cca", dims = 1:10) %>%
  FindClusters(resolution = 0.5, cluster.name = "integrated.clusters") %>%
  RunUMAP(dims = 1:10, reduction = "integrated.cca", reduction.name = "integrated.umap")

# Step 6: Save the integrated object for further analysis
saveRDS(ser_int, file = "data/analysis/integrated_sn_human_PVAT_seurat_object.rds")

# Visualisation
print(ser_int[['pca']], dims = 1:5, nfeatures = 5)
VizDimLoadings(ser_int, dims = 1:2, reduction = 'pca')
DimPlot(ser_int, reduction = 'pca')
DimHeatmap(ser_int, dims = 1:15, cells = 500, balanced = TRUE)
ElbowPlot(ser_int)
head(Idents(ser_int), 5)
DimPlot(ser_int, reduction = 'umap')

# Visualise umap coloured on sample origin
ser_int$sample_origin <- sub("_.*", "", colnames(ser_int))
Idents(ser_int) <- "sample_origin"
DimPlot(ser_int, reduction = "umap", group.by = "sample_origin", label = TRUE, repel = TRUE) +
  ggtitle("UMAP Colored by Sample Origin") +
  theme_minimal()
# Visualise PCA coloured by sample origin
DimPlot(ser_int, reduction = "pca", group.by = "sample_origin", label = TRUE)

##### Additional code ##### 
# Add sample information to each Seurat object
#human1$sample <- "human1"
#human2$sample <- "human2"
#human3$sample <- "human3"

# Merge the three Seurat objects
#human_int <- merge(human1, y = c(human2, human3), add.cell.ids = c("human1", "human2", "human3"))

# Split the RNA assay into layers by sample
#human_int[["RNA"]] <- split(human_int[["RNA"]], f = human_int$sample)

# Integrate layers
#human_int <- IntegrateLayers(
#  object = human_int,
#  method = HarmonyIntegration, 
#  orig.reduction = "pca", 
#  new.reduction = "harmony",
#  verbose = FALSE
#)
