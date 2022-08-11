# if (!require("BiocManager", quietly = TRUE))
#   install.packages("BiocManager")
# 
# BiocManager::install("SingleCellExperiment")
# BiocManager::install("batchelor")
# install.packages("harmony")
# install.packages("Seurat")

# if(!require('cowplot')) {
#   install.packages('cowplot')
# }
# library(devtools)
# install_github("immunogenomics/harmony")

library(SingleCellExperiment)
library(harmony)
library(Seurat)
library('cowplot')

#Existing dataframes after gene expression extraction are:
#counts_Xin and counts_Lawlor 
counts_Xin_Lawlor = cbind(counts_Xin,counts_Lawlor)

batch_Xin = rep("Xin", dim(counts_Xin)[1])
batch_Lawlor = rep("Lawlor", dim(counts_Lawlor)[1])
batch_Xin_Lawlor = data.frame(cbind(batch_Xin,batch_Lawlor))

#Using Seurat
Sobj <- CreateSeuratObject(counts = counts_Xin_Lawlor, project = "Xin_Lawlor", min.cells = 5) %>%
  Seurat::NormalizeData(verbose = FALSE) %>%
  FindVariableFeatures(selection.method = "vst", nfeatures = 2000) %>% 
  ScaleData(verbose = FALSE) %>% 
  RunPCA(npcs = 20, verbose = FALSE)

Sobj@meta.data$batch <- c(rep("Xin", dim(counts_Xin)[2]), rep("Lawlor", dim(counts_Lawlor)[2]))

options(repr.plot.height = 5, repr.plot.width = 12)
p1 <- DimPlot(object = Sobj, reduction = "pca", pt.size = .1, group.by = "batch")

options(repr.plot.height = 2.5, repr.plot.width = 6)
Sobj <- Sobj %>% 
  RunHarmony("batch", plot_convergence = TRUE)

options(repr.plot.height = 5, repr.plot.width = 12)
p2 <- DimPlot(object = Sobj, reduction = "harmony", pt.size = .1, group.by = "batch")

plot_grid(p1,p2)

Sobj <- Sobj %>% 
  RunUMAP(reduction = "harmony", dims = 1:20) %>% 
  FindNeighbors(reduction = "harmony", dims = 1:20) %>% 
  FindClusters(resolution = 0.5) %>% 
  identity()

options(repr.plot.height = 4, repr.plot.width = 10)
DimPlot(Sobj, reduction = "umap", group.by = "batch", pt.size = .1, split.by = 'batch')

options(repr.plot.height = 4, repr.plot.width = 6)
p4 <- DimPlot(Sobj, reduction = "umap", group.by = "batch",pt.size = .1)

options(repr.plot.height = 4, repr.plot.width = 6)
p3 <- DimPlot(Sobj, reduction = "umap", label = TRUE, pt.size = .1)

plot_grid(p3,p4)

Sobj.harmony.embed <- Embeddings(Sobj,"harmony")