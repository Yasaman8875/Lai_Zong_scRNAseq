library(Seurat)
library("tidyverse")
library("clustree")
library("future")
library("unixtools")
library(ggplot2)
library(patchwork)
library("scProportionTest")

options(future.globals.maxSize = 10000 * 1024 ^2)
plan("multiprocess", workers = 2)

samples_10x <- list(DMSO = file.path("DMSO", "outs", "filtered_feature_bc_matrix"),
                    GSK126 = file.path("EZH2i", "outs", "filtered_feature_bc_matrix"),
                    NSC23766 = file.path("RACi", "outs", "filtered_feature_bc_matrix"),
                    Combo = file.path("EZH2-RACi", "outs", "filtered_feature_bc_matrix")
)
counts_10X <- map(samples_10x, Read10X)
seurat_obj <- imap(counts_10X, ~ CreateSeuratObject(
  counts = .x, project = .y, min.cells = 10, min.features = 250
))

## Cell Quality Control
## ----------

## Add mitochondrial percentage.

seurat_obj <- map(seurat_obj, function(x) {
  x[["percent.mt"]] <- PercentageFeatureSet(x, "^MT-")
  return(x)
})

## Add complexity scores.

seurat_obj <- map(seurat_obj, function(x) {
  x$log10GenesPerUMI <- log10(x$nFeature_RNA) / log10(x$nCount_RNA)
  return(x)
})

## Prepare plots

mt_data <- map(seurat_obj, function(x) {
  meta_data <- x@meta.data
  meta_data <- as_tibble(meta_data, .name_repair = "unique")
  meta_data$cells <- rownames(meta_data) 
  meta_data <- meta_data %>% rename(nUmi = nCount_RNA, nGene = nFeature_RNA)
  return(meta_data)
})

mt_data$DMSO$sample <- "DMSO"
mt_data$GSK126$sample <- "GSK126"
mt_data$NSC23766$sample <- "NSC23766"
mt_data$Combo$sample <- "Combo"

mt_data <- bind_rows(mt_data)

p <- ggplot(mt_data, aes(x = percent.mt, y = nGene)) +
  geom_point(size = 0.1) +
  facet_wrap(. ~ orig.ident, ncol = 3, scales = "free")

if (!dir.exists(file.path("results", "preprocessing"))) {
  dir.create(file.path("results", "preprocessing"))
}

pdf(file.path("results", "preprocessing", "mt_content.pdf"), height = 10, width = 10)
p
dev.off()

p <- ggplot(mt_data, aes(x=sample, fill=sample)) + 
  geom_bar() + theme_classic() +
  theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust=1)) +
  theme(plot.title = element_text(hjust=0.5, face="bold")) +
  ggtitle("NCells")

pdf(file.path("results", "preprocessing", "num_cells.pdf"))
p
dev.off()

p <- ggplot(mt_data, aes(color=sample, x=nUmi, fill=sample)) +
  geom_density(alpha = 0.2) +
  scale_x_log10() +
  theme_classic() +
  ylab("Cell density") +
  geom_vline(xintercept = 500)
pdf(file.path("results", "preprocessing", "umi_counts.pdf"), height = 10, width = 10)
p
dev.off()

p <- ggplot(mt_data, aes(color=sample, x = nGene, fill= sample))  +
  geom_density(alpha = 0.2) +
  theme_classic() +
  scale_x_log10() +
  geom_vline(xintercept = 300)
pdf(file.path("results", "preprocessing", "genes_per_cell.pdf"), height = 10, width = 10)
p
dev.off()


p <- ggplot(mt_data, aes(x = sample, y = log10(nGene), fill=sample))+
  geom_boxplot() +
  theme_classic() +
  theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust=1)) +
  theme(plot.title = element_text(hjust=0.5, face="bold")) +
  ggtitle("NCells vs NGenes")

pdf(file.path("results", "preprocessing", "nGenes.pdf"), height = 10, width = 10)
p
dev.off()

p <- ggplot(mt_data, aes(x = nUmi, y = nGene, color = percent.mt)) +
  geom_point() +
  scale_colour_gradient(low = "gray90", high = "black") +
  stat_smooth(method=lm) +
  scale_x_log10() +
  scale_y_log10() +
  theme_classic() +
  geom_vline(xintercept = 500) +
  geom_hline(yintercept = 250) +
  facet_wrap(~sample)
pdf(file.path("results", "preprocessing", "umi_v_genes.pdf"), height = 10, width = 10)
p
dev.off()

p <- ggplot(mt_data, aes(color=sample, x = percent.mt, fill=sample)) +
  geom_density(alpha = 0.2) +
  scale_x_log10() +
  theme_classic() +
  geom_vline(xintercept = 20)
pdf(file.path("results", "preprocessing", "mt_density.pdf"), height = 10, width = 10)
p
dev.off()

p <- ggplot(mt_data, aes(x=log10GenesPerUMI, color = sample, fill=sample)) +
  geom_density(alpha = 0.2) +
  theme_classic() +
  geom_vline(xintercept = 0.8)

pdf(file.path("results", "preprocessing", "complexity.pdf"), height = 10, width = 10)
p
dev.off()

seurat_obj <- map(seurat_obj, function(x) {
  x <- subset(x, subset = percent.mt <= 20
              & nFeature_RNA > 1000 & log10GenesPerUMI > 0.8)
}
)
seurat_obj <- map(seurat_obj, function(x) {
  x <- CellCycleScoring(
    x, s.features = cc.genes$s.genes, g2m.features = cc.genes$g2m.genes,
    set.ident = TRUE
  )
  return(x)
})
seurat_obj <- map(seurat_obj, function(x) {
  x$CC.Difference <- x$S.Score - x$G2M.Score
  return(x)
})
seurat_obj <- map(seurat_obj, function(x) {
  x <- SCTransform(x, vars.to.regress = c("CC.Difference", "percent.mt"), verbose = FALSE)
  return(x)
})

saveRDS(seurat_obj, file.path("results", "r_objects", "seurat_obj.RDS"))
integration_features <- SelectIntegrationFeatures(seurat_obj, nfeatures = 3000)
seurat_obj <- PrepSCTIntegration(seurat_obj, anchor.features = integration_features)

integration_anchors <- FindIntegrationAnchors(
  seurat_obj, normalization.method = "SCT", anchor.features = integration_features)
seurat_integrated <- IntegrateData(integration_anchors, normalization.method = "SCT")

seurat_integrated <- RunPCA(seurat_integrated, npcs = 100)
saveRDS(seurat_integrated, file.path("results", "r_objects", "seurat_integrated.RDS"))
p <- ElbowPlot(seurat_integrated, ndims = 100)

if (!dir.exists(file.path("results", "clustering"))) {
  dir.create(file.path("results", "clustering"))
}

pdf(file.path("results", "clustering", "pca_elbow_plot.pdf"), height = 5, width = 5)
p
dev.off()
seurat_integrated <- JackStraw(seurat_integrated, num.replicate = 100, dims = 50)
seurat_integrated <- ScoreJackStraw(seurat_integrated, dims = 1:45)
p <- JackStrawPlot(seurat_integrated, dims = 1:45)
pdf(file.path("results", "clustering", "JackStrawPlot.pdf"), height = 10, width = 12)
p
dev.off()

seurat_integrated <- FindNeighbors(seurat_integrated, dims = 1:40)
seurat_integrated <- FindClusters(
  seurat_integrated, resolution = seq(0.2, 1.2, 0.1),
  algorithm = 4, method = "igraph", weights = TRUE
)
## Plotting a cluster tree.

p <- clustree(seurat_integrated, prefix = "integrated_snn_res.")

pdf(file.path("results", "clustering", "cluster_tree.pdf"), height = 16, width = 12)
p
dev.off()

if (!dir.exists("tempdir")) dir.create("tempdir")
set.tempdir("tempdir")

Idents(seurat_integrated) <- "integrated_snn_res.0.3"
seurat_integrated <- RunUMAP(seurat_integrated, dims = 1:40)

seurat_integrated$orig.ident <- factor(
  seurat_integrated$orig.ident, levels = c("DMSO", "GSK126", "NSC23766", "Combo"))

p <- DimPlot(seurat_integrated, group.by = "ident", split.by = "orig.ident", ncol = 2)

pdf(file.path("results", "clustering", "clusters.pdf"), height = 12, width = 10)
p
dev.off()

# Integrated clusters
p <- DimPlot(seurat_integrated)
pdf(file.path("results", "clustering", "integrated_clusters.pdf"), height = 10, width = 10)
p
dev.off()

# QC Clusters
metrics <-  c("nCount_RNA", "nFeature_RNA", "S.Score", "G2M.Score", "percent.mt")


p <- FeaturePlot(seurat_integrated, 
                 reduction = "umap", 
                 features = metrics,
                 pt.size = 0.4, 
                 sort.cell = TRUE,
                 min.cutoff = 'q10',
                 label = TRUE)
pdf(file.path("results", "clustering", "QCFeatures.pdf"))
p
dev.off()
## Save Integrated Data
## ----------

saveRDS(seurat_integrated, file.path("results", "r_objects", "seurat_integrated.RDS"))
