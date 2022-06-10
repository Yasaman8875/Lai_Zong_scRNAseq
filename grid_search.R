library(Seurat)
library("tidyverse")
library("clustree")
library("future")
library("unixtools")
library(ggplot2)
library(patchwork)
library("scProportionTest")

options(future.globals.maxSize = 10000 * 1024 ^2)
plan("multiprocess", workers = 6)

seurat_obj <- readRDS(file.path("results", "r_objects", "seurat_obj.RDS"))
integration_features <- SelectIntegrationFeatures(seurat_obj, nfeatures = 3000)
seurat_obj <- PrepSCTIntegration(seurat_obj, anchor.features = integration_features)

integration_anchors <- FindIntegrationAnchors(
    seurat_obj, normalization.method = "SCT", anchor.features = integration_features)
seurat_integrated <- IntegrateData(integration_anchors, normalization.method = "SCT")

## Original analysis used 40 PCA dimensions, k.param = 20, prune.SNN = 1/15, res = 0.3

if (!dir.exists(file.path("results", "grid"))) {
    dir.create(file.path("results", "grid"))
}

seurat_integrated <- RunPCA(seurat_integrated, npcs = 50)
saveRDS(seurat_integrated, file.path("results", "grid", "seurat_integrated.RDS"))


# https://rpubs.com/kinsimon96/mpeg_analysis


seurat_integrated <- FindNeighbors(seurat_integrated, dims = 1:40, k.param = 15)
seurat_integrated <- FindClusters(
    seurat_integrated, resolution = seq(0.2, 1.2, 0.1),
    algorithm = 4, method = "igraph", weights = TRUE
)
saveRDS(seurat_integrated, file.path("results", "grid", "seurat_integrated_k15.RDS"))
## Plotting a cluster tree.

p <- clustree(seurat_integrated, prefix = "integrated_snn_res.")

pdf(file.path("results", "grid", "clustree_k15.pdf"), height = 16, width = 12)
p
dev.off()

# Reload object to calculate another iteration
seurat_integrated <- readRDS(file.path("results", "grid", "seurat_integrated.RDS"))
seurat_integrated <- FindNeighbors(seurat_integrated, dims = 1:40, k.param = 25)
seurat_integrated <- FindClusters(
    seurat_integrated, resolution = seq(0.2, 1.2, 0.1),
    algorithm = 4, method = "igraph", weights = TRUE
)
saveRDS(seurat_integrated, file.path("results", "grid", "seurat_integrated_k25.RDS"))
## Plotting a cluster tree.

p <- clustree(seurat_integrated, prefix = "integrated_snn_res.")

pdf(file.path("results", "grid", "clustree_k25.pdf"), height = 16, width = 12)
p
dev.off()

# Reload object to calculate another iteration
seurat_integrated <- readRDS(file.path("results", "grid", "seurat_integrated.RDS"))
seurat_integrated <- FindNeighbors(seurat_integrated, dims = 1:40, k.param = 25, prune.SNN = 1/20)
seurat_integrated <- FindClusters(
    seurat_integrated, resolution = seq(0.2, 1.2, 0.1),
    algorithm = 4, method = "igraph", weights = TRUE
)
saveRDS(seurat_integrated, file.path("results", "grid", "seurat_integrated_k25p20.RDS"))
## Plotting a cluster tree.

p <- clustree(seurat_integrated, prefix = "integrated_snn_res.")

pdf(file.path("results", "grid", "clustree_k25_p20.pdf"), height = 16, width = 12)
p
dev.off()

# Reload object to calculate another iteration
seurat_integrated <- readRDS(file.path("results", "grid", "seurat_integrated.RDS"))
seurat_integrated <- FindNeighbors(seurat_integrated, dims = 1:40, k.param = 25, prune.SNN = 2/15)
seurat_integrated <- FindClusters(
    seurat_integrated, resolution = seq(0.2, 1.2, 0.1),
    algorithm = 4, method = "igraph", weights = TRUE
)
saveRDS(seurat_integrated, file.path("results", "grid", "seurat_integrated_k25p2.RDS"))
## Plotting a cluster tree.

p <- clustree(seurat_integrated, prefix = "integrated_snn_res.")

pdf(file.path("results", "grid", "clustree_k25_p2.pdf"), height = 16, width = 12)
p
dev.off()

# Reload object to calculate another iteration
seurat_integrated <- readRDS(file.path("results", "grid", "seurat_integrated.RDS"))
seurat_integrated <- FindNeighbors(seurat_integrated, dims = 1:40, k.param = 15, prune.SNN = 1/20)
seurat_integrated <- FindClusters(
    seurat_integrated, resolution = seq(0.2, 1.2, 0.1),
    algorithm = 4, method = "igraph", weights = TRUE
)
saveRDS(seurat_integrated, file.path("results", "grid", "seurat_integrated_k15p20.RDS"))
## Plotting a cluster tree.

p <- clustree(seurat_integrated, prefix = "integrated_snn_res.")

pdf(file.path("results", "grid", "clustree_k15_p20.pdf"), height = 16, width = 12)
p
dev.off()

# Reload object to calculate another iteration
seurat_integrated <- readRDS(file.path("results", "grid", "seurat_integrated.RDS"))
seurat_integrated <- FindNeighbors(seurat_integrated, dims = 1:40, k.param = 15, prune.SNN = 2/15)
seurat_integrated <- FindClusters(
    seurat_integrated, resolution = seq(0.2, 1.2, 0.1),
    algorithm = 4, method = "igraph", weights = TRUE
)
saveRDS(seurat_integrated, file.path("results", "grid", "seurat_integrated_k15p2.RDS"))
## Plotting a cluster tree.

p <- clustree(seurat_integrated, prefix = "integrated_snn_res.")

pdf(file.path("results", "grid", "clustree_k15_p2.pdf"), height = 16, width = 12)
p
dev.off()
