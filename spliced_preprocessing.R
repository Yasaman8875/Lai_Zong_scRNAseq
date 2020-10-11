library("Seurat")
library("tidyverse")
library("clustree")
library("future")
library("unixtools")
library("loomR")
library("velocyto.R")
library("SeuratWrappers")

options(future.globals.maxSize = 10000 * 1024 ^2)
plan("multiprocess", workers = 4)

velocyto_samples <- list(
  DMSO = file.path("DMSO", "velocyto", "DMSO.loom"),
  EZH2i = file.path("EZH2i", "velocyto", "EZH2i.loom"),
  RACi = file.path("RACi", "velocyto", "RACi.loom"),
  Combo = file.path("EZH2-RACi", "velocyto", "EZH2-RACi.loom")
)

## Read in data.

velocyto_counts <- map(velocyto_samples, ReadVelocity)

## Convert to seurat object.

seurat_obj <- map(velocyto_counts, as.Seurat)

seurat_obj <- imap(seurat_obj, function(x, y) {
  x@meta.data$orig.ident <- y
  return(x)
})

## Add complexity scores.

seurat_obj <- map(seurat_obj, function(x) {
  x$log10GenesPerUMI_spliced <- log10(x$nFeature_spliced) / log10(x$nCount_spliced)
  return(x)
})

mt_data <- map(seurat_obj, function(x) {
  meta_data <- x@meta.data
  meta_data <- as_tibble(meta_data, .name_repair = "unique")
  meta_data$cells <- rownames(meta_data) 
  meta_data <- meta_data %>% rename(nUmi = nCount_spliced, nGene = nFeature_spliced)
  return(meta_data)
})

mt_data$DMSO$sample <- "DMSO"
mt_data$EZH2$sample <- "EZH2"
mt_data$RACi$sample <- "RACi"
mt_data$Combo$sample <- "Combo"

mt_data <- bind_rows(mt_data)

p <- ggplot(mt_data, aes(color=sample, x=nUmi, fill=sample)) +
  geom_density(alpha = 0.2) +
  scale_x_log10() +
  theme_classic() +
  ylab("Cell density") +
  geom_vline(xintercept = 500)
pdf(file.path("results", "preprocessing", "umi_counts_spliced.pdf"), height = 10, width = 10)
p
dev.off()

p <- ggplot(mt_data, aes(color=sample, x = nGene, fill= sample))  +
  geom_density(alpha = 0.2) +
  theme_classic() +
  scale_x_log10() +
  geom_vline(xintercept = 300)
pdf(file.path("results", "preprocessing", "spliced_genes_per_cell.pdf"), height = 10, width = 10)
p
dev.off()

p <- ggplot(mt_data, aes(x = sample, y = log10(nGene), fill=sample))+
  geom_boxplot() +
  theme_classic() +
  theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust=1)) +
  theme(plot.title = element_text(hjust=0.5, face="bold")) +
  ggtitle("NCells vs NGenes")

pdf(file.path("results", "preprocessing", "nGenes_spliced.pdf"), height = 10, width = 10)
p
dev.off()

p <- ggplot(mt_data, aes(x=log10GenesPerUMI_spliced, color = sample, fill=sample)) +
  geom_density(alpha = 0.2) +
  theme_classic() +
  geom_vline(xintercept = 0.8)

pdf(file.path("results", "preprocessing", "complexity_spliced.pdf"), height = 10, width = 10)
p
dev.off()

saveRDS(seurat_obj, file.path("results", "r_objects", "seurat_obj_spliced.RDS"))

## Filter data
seurat_obj <- map(seurat_obj, function(x) {
  x <- subset(x, subset = log10GenesPerUMI_spliced > 0.8 
              & nFeature_spliced > 200 & nCount_spliced > 350)
}
)

## SCTransform to normalize data.

seurat_obj <- map(seurat_obj, ~ SCTransform(., assay = "spliced"))

saveRDS(seurat_obj, file.path("results", "r_objects", "seurat_obj_spliced.RDS"))
## Prepare for integration.

integration_features <- SelectIntegrationFeatures(seurat_obj, nfeatures = 3000)
seurat_obj <- PrepSCTIntegration(seurat_obj, anchor.features = integration_features)
integration_anchors <- FindIntegrationAnchors(
  seurat_obj, normalization.method = "SCT", anchor.features = integration_features)
seurat_integrated <- IntegrateData(integration_anchors, normalization.method = "SCT")

# Dimension Reduction
seurat_integrated <- RunPCA(seurat_integrated, npcs = 50)
p <- ElbowPlot(seurat_integrated, ndims = 50)

if (!dir.exists(file.path("results", "clustering"))) {
  dir.create(file.path("results", "clustering"))
}

pdf(file.path("results", "clustering", "pca_elbow_plot_spliced.pdf"), height = 5, width = 5)
p
dev.off()

seurat_integrated <- JackStraw(seurat_integrated, num.replicate = 100, dims = 30)
seurat_integrated <- ScoreJackStraw(seurat_integrated, dims = 1:30)
p <- JackStrawPlot(seurat_integrated, dims = 1:30)
pdf(file.path("results", "clustering", "JackStrawPlot_spliced.pdf"), height = 10, width = 12)
p
dev.off()

seurat_integrated <- FindNeighbors(seurat_integrated, dims = 1:25)
seurat_integrated <- FindClusters(
  seurat_integrated, resolution = seq(0.2, 1.2, 0.1),
  method = "igraph", algorithm = 4, weights = TRUE
)

saveRDS(seurat_integrated, file.path("results", "r_objects", "seurat_integrated_spliced.RDS"))

## Plotting a cluster tree.

p <- clustree(seurat_integrated, prefix = "integrated_snn_res.")

pdf(file.path("results", "clustering", "cluster_tree_spliced.pdf"), height = 16, width = 12)
p
dev.off()
## Switch identity to a presumptive good clustering resolution.

Idents(seurat_integrated) <- "integrated_snn_res.0.3"

## UMAP dimension reduction for visualization.

if (!dir.exists("tempdir")) dir.create("tempdir")
set.tempdir("tempdir")

seurat_integrated <- RunUMAP(seurat_integrated, dims = 1:30)

saveRDS(seurat_integrated, file.path("results", "r_objects", "seurat_integrated_spliced.RDS"))

## Plot Clusters.

p <- DimPlot(seurat_integrated, group.by = "ident", split.by = "orig.ident", ncol = 2)

pdf(file.path("results", "clustering", "clusters_spliced.pdf"), height = 12, width = 10)
p
dev.off()

