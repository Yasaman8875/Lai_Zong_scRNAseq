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
                    EZH2i = file.path("EZH2i", "outs", "filtered_feature_bc_matrix"),
                    RACi = file.path("RACi", "outs", "filtered_feature_bc_matrix"),
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
  x <- SCTransform(x, vars.to.regress = "CC.Difference", verbose = FALSE)
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

Idents(seurat_integrated) <- "integrated_snn_res.0.7"
seurat_integrated <- RunUMAP(seurat_integrated, dims = 1:40)
p <- DimPlot(seurat_integrated, group.by = "ident", split.by = "orig.ident", ncol = 2)

pdf(file.path("results", "clustering", "clusters.pdf"), height = 12, width = 10)
p
dev.off()

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

## Marker Analysis

if (!dir.exists(file.path("results", "markers"))) {
  dir.create(file.path("results", "markers"))
}

get_conserved <- function(cluster){
  FindConservedMarkers(seurat_integrated,
                       ident.1 = cluster,
                       grouping.var = "orig.ident",
                       only.pos = TRUE, 
                       assay = "SCT", 
                       slot = "data") %>%
    rownames_to_column(var = "gene") %>%
    left_join(y = unique(annotations[, c("gene_name", "description")]),
              by = c("gene" = "gene_name")) %>%
    cbind(cluster_id = cluster, .)
}
annotations <- read.csv("annotation.csv")

conserved_markers <- map_dfr(1:20, get_conserved)

top10 <- conserved_markers %>% 
  mutate(avg_fc = (control_avg_logFC + crispr_avg_logFC) /2) %>% 
  group_by(cluster_id) %>% 
  top_n(n = 10, 
        wt = avg_fc)
write.table(top10, file.path("results", "markers", "top10.csv"), sep=",", col.names=TRUE, row.names=FALSE, quote=FALSE)
