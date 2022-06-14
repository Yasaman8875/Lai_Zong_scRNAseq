library(Seurat)
library("tidyverse")
library("clustree")
library("future")
library("unixtools")
library(ggplot2)
library(patchwork)
library("scProportionTest")
library(factoextra)

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


seurat_integrated <- readRDS(file.path("results", "grid", "seurat_integrated.RDS"))
seurat_integrated <- FindNeighbors(seurat_integrated, dims = 1:40, k.param = 18)
seurat_integrated <- FindClusters(
    seurat_integrated, resolution = seq(0.3, 0.4, 0.1),
    algorithm = 4, method = "igraph", weights = TRUE
)
saveRDS(seurat_integrated, file.path("results", "grid", "seurat_integrated_k18.RDS"))

seurat_integrated <- readRDS(file.path("results", "grid", "seurat_integrated.RDS"))
seurat_integrated <- FindNeighbors(seurat_integrated, dims = 1:40, k.param = 18, prune.SNN = 1/12)
seurat_integrated <- FindClusters(
    seurat_integrated, resolution = seq(0.3, 0.4, 0.1),
    algorithm = 4, method = "igraph", weights = TRUE
)
saveRDS(seurat_integrated, file.path("results", "grid", "seurat_integrated_k18p12.RDS"))

seurat_integrated <- readRDS(file.path("results", "grid", "seurat_integrated.RDS"))
seurat_integrated <- FindNeighbors(seurat_integrated, dims = 1:40, k.param = 18, prune.SNN = 1/17)
seurat_integrated <- FindClusters(
    seurat_integrated, resolution = seq(0.3, 0.4, 0.1),
    algorithm = 4, method = "igraph", weights = TRUE
)
saveRDS(seurat_integrated, file.path("results", "grid", "seurat_integrated_k18p17.RDS"))

seurat_integrated <- readRDS(file.path("results", "grid", "seurat_integrated.RDS"))
seurat_integrated <- FindNeighbors(seurat_integrated, dims = 1:40, k.param = 22)
seurat_integrated <- FindClusters(
    seurat_integrated, resolution = seq(0.3, 0.4, 0.1),
    algorithm = 4, method = "igraph", weights = TRUE
)
saveRDS(seurat_integrated, file.path("results", "grid", "seurat_integrated_k22.RDS"))

seurat_integrated <- readRDS(file.path("results", "grid", "seurat_integrated.RDS"))
seurat_integrated <- FindNeighbors(seurat_integrated, dims = 1:40, k.param = 22, prune.SNN = 1/12)
seurat_integrated <- FindClusters(
    seurat_integrated, resolution = seq(0.3, 0.4, 0.1),
    algorithm = 4, method = "igraph", weights = TRUE
)
saveRDS(seurat_integrated, file.path("results", "grid", "seurat_integrated_k22p12.RDS"))

seurat_integrated <- readRDS(file.path("results", "grid", "seurat_integrated.RDS"))
seurat_integrated <- FindNeighbors(seurat_integrated, dims = 1:40, k.param = 22, prune.SNN = 1/17)
seurat_integrated <- FindClusters(
    seurat_integrated, resolution = seq(0.3, 0.4, 0.1),
    algorithm = 4, method = "igraph", weights = TRUE
)
saveRDS(seurat_integrated, file.path("results", "grid", "seurat_integrated_k22p17.RDS"))

seurat_integrated <- readRDS(file.path("results", "grid", "seurat_integrated.RDS"))
seurat_integrated <- FindNeighbors(seurat_integrated, dims = 1:40, prune.SNN = 1/12)
seurat_integrated <- FindClusters(
    seurat_integrated, resolution = seq(0.3, 0.4, 0.1),
    algorithm = 4, method = "igraph", weights = TRUE
)
saveRDS(seurat_integrated, file.path("results", "grid", "seurat_integrated_p12.RDS"))

seurat_integrated <- readRDS(file.path("results", "grid", "seurat_integrated.RDS"))
seurat_integrated <- FindNeighbors(seurat_integrated, dims = 1:40, prune.SNN = 1/17)
seurat_integrated <- FindClusters(
    seurat_integrated, resolution = seq(0.3, 0.4, 0.1),
    algorithm = 4, method = "igraph", weights = TRUE
)
saveRDS(seurat_integrated, file.path("results", "grid", "seurat_integrated_p17.RDS"))

silhouette_plot <- function(seurat_integrated) {
    # https://rpubs.com/kinsimon96/mpeg_analysis
    #Compute distance matrix to UMAP coordinates
    distance_matrix <- dist(Embeddings(seurat_integrated[['umap']])[, 1:2])
    
    #Isolate cluster nameq
    clusters <- seurat_integrated@active.ident
    
    #Compute silhouette score for each cluster
    silhouette <- silhouette(as.numeric(clusters), dist = distance_matrix)
    seurat_integrated@meta.data$silhouette_score <- silhouette[,3]
    
    mean_silhouette_score <- mean(seurat_integrated@meta.data$silhouette_score)
    p <- seurat_integrated@meta.data %>%
        mutate(barcode = rownames(.)) %>%
        arrange(seurat_clusters,-silhouette_score) %>%
        mutate(barcode = factor(barcode, levels = barcode)) %>%
        ggplot() +
        geom_col(aes(barcode, silhouette_score, fill = seurat_clusters), show.legend = FALSE) +
        geom_hline(yintercept = mean_silhouette_score, color = 'red', linetype = 'dashed') +
        scale_x_discrete(name = 'Cells') +
        scale_y_continuous(name = 'Silhouette score') +
        theme_bw() +
        theme(
            axis.title.x = element_blank(),
            axis.text.x = element_blank(),
            axis.ticks.x = element_blank(),
            panel.grid.major = element_blank(),
            panel.grid.minor = element_blank()
        )
    return(p)
}

seurat_integrated <- readRDS(file.path("results", "r_objects", "seurat_integrated.RDS"))
Idents(seurat_integrated) <- "integrated_snn_res.0.3"
DefaultAssay(seurat_integrated) <- "SCT"

p <- silhouette_score(seurat_integrated)
pdf(file.path("results", "grid", "base_sil.pdf"), height = 12, width = 12)
p
dev.off()

# Get a more detailed silhouette plot 

distance_matrix <- dist(Embeddings(seurat_integrated[['umap']])[, 1:2])
clusters <- seurat_integrated@active.ident
silhouette <- silhouette(as.numeric(clusters), dist = distance_matrix)
p <- fviz_silhouette(silhouette)
png(file.path("results", "grid", "base_sil.png"))

seurat_integrated <- readRDS(file.path("results", "grid", "seurat_integrated_k15.RDS"))
Idents(seurat_integrated) <- "integrated_snn_res.0.3"
seurat_integrated <- RunUMAP(seurat_integrated, dims = 1:40)

p <- silhouette_score(seurat_integrated)
pdf(file.path("results", "grid", "k15r3_sil.pdf"), height = 12, width = 12)
p
dev.off()

Idents(seurat_integrated) <- "integrated_snn_res.0.4"
seurat_integrated <- RunUMAP(seurat_integrated, dims = 1:40)

p <- silhouette_score(seurat_integrated)
pdf(file.path("results", "grid", "k15r4_sil.pdf"), height = 12, width = 12)
p
dev.off()

seurat_integrated <- readRDS(file.path("results", "grid", "seurat_integrated_k15p20.RDS"))
Idents(seurat_integrated) <- "integrated_snn_res.0.3"
seurat_integrated <- RunUMAP(seurat_integrated, dims = 1:40)

p <- silhouette_score(seurat_integrated)
pdf(file.path("results", "grid", "k15p20r3_sil.pdf"), height = 12, width = 12)
p
dev.off()

Idents(seurat_integrated) <- "integrated_snn_res.0.4"
seurat_integrated <- RunUMAP(seurat_integrated, dims = 1:40)

p <- silhouette_score(seurat_integrated)
pdf(file.path("results", "grid", "k15p20r4_sil.pdf"), height = 12, width = 12)
p
dev.off()

seurat_integrated <- readRDS(file.path("results", "grid", "seurat_integrated_k15p2.RDS"))
Idents(seurat_integrated) <- "integrated_snn_res.0.3"
seurat_integrated <- RunUMAP(seurat_integrated, dims = 1:40)

p <- silhouette_score(seurat_integrated)
pdf(file.path("results", "grid", "k15p2r3_sil.pdf"), height = 12, width = 12)
p
dev.off()

Idents(seurat_integrated) <- "integrated_snn_res.0.4"
seurat_integrated <- RunUMAP(seurat_integrated, dims = 1:40)

p <- silhouette_score(seurat_integrated)
pdf(file.path("results", "grid", "k15p2r4_sil.pdf"), height = 12, width = 12)
p
dev.off()

seurat_integrated <- readRDS(file.path("results", "grid", "seurat_integrated_k25.RDS"))
Idents(seurat_integrated) <- "integrated_snn_res.0.3"
seurat_integrated <- RunUMAP(seurat_integrated, dims = 1:40)

p <- silhouette_score(seurat_integrated)
pdf(file.path("results", "grid", "k25r3_sil.pdf"), height = 12, width = 12)
p
dev.off()

Idents(seurat_integrated) <- "integrated_snn_res.0.4"
seurat_integrated <- RunUMAP(seurat_integrated, dims = 1:40)

p <- silhouette_score(seurat_integrated)
pdf(file.path("results", "grid", "k25r4_sil.pdf"), height = 12, width = 12)
p
dev.off()

seurat_integrated <- readRDS(file.path("results", "grid", "seurat_integrated_k25p20.RDS"))
Idents(seurat_integrated) <- "integrated_snn_res.0.3"
seurat_integrated <- RunUMAP(seurat_integrated, dims = 1:40)

p <- silhouette_score(seurat_integrated)
pdf(file.path("results", "grid", "k25p20r3_sil.pdf"), height = 12, width = 12)
p
dev.off()

Idents(seurat_integrated) <- "integrated_snn_res.0.4"
seurat_integrated <- RunUMAP(seurat_integrated, dims = 1:40)

p <- silhouette_score(seurat_integrated)
pdf(file.path("results", "grid", "k25p20r4_sil.pdf"), height = 12, width = 12)
p
dev.off()

seurat_integrated <- readRDS(file.path("results", "grid", "seurat_integrated_k25p2.RDS"))
Idents(seurat_integrated) <- "integrated_snn_res.0.3"
seurat_integrated <- RunUMAP(seurat_integrated, dims = 1:40)

p <- silhouette_score(seurat_integrated)
pdf(file.path("results", "grid", "k25p2r3_sil.pdf"), height = 12, width = 12)
p
dev.off()

Idents(seurat_integrated) <- "integrated_snn_res.0.4"
seurat_integrated <- RunUMAP(seurat_integrated, dims = 1:40)

p <- silhouette_score(seurat_integrated)
pdf(file.path("results", "grid", "k25p2r4_sil.pdf"), height = 12, width = 12)
p
dev.off()


seurat_integrated <- readRDS(file.path("results", "grid", "seurat_integrated_p17.RDS"))
Idents(seurat_integrated) <- "integrated_snn_res.0.3"
seurat_integrated <- RunUMAP(seurat_integrated, dims = 1:40)

p <- silhouette_score(seurat_integrated)
pdf(file.path("results", "grid", "p17r3_sil.pdf"), height = 12, width = 12)
p
dev.off()

Idents(seurat_integrated) <- "integrated_snn_res.0.4"
seurat_integrated <- RunUMAP(seurat_integrated, dims = 1:40)

p <- silhouette_score(seurat_integrated)
pdf(file.path("results", "grid", "p17r4_sil.pdf"), height = 12, width = 12)
p
dev.off()

seurat_integrated <- readRDS(file.path("results", "grid", "seurat_integrated_k18.RDS"))
Idents(seurat_integrated) <- "integrated_snn_res.0.3"
seurat_integrated <- RunUMAP(seurat_integrated, dims = 1:40)

p <- silhouette_score(seurat_integrated)
pdf(file.path("results", "grid", "k18r3_sil.pdf"), height = 12, width = 12)
p
dev.off()

Idents(seurat_integrated) <- "integrated_snn_res.0.4"
seurat_integrated <- RunUMAP(seurat_integrated, dims = 1:40)

p <- silhouette_score(seurat_integrated)
pdf(file.path("results", "grid", "k18r4_sil.pdf"), height = 12, width = 12)
p
dev.off()

seurat_integrated <- readRDS(file.path("results", "grid", "seurat_integrated_k18p12.RDS"))
Idents(seurat_integrated) <- "integrated_snn_res.0.3"
seurat_integrated <- RunUMAP(seurat_integrated, dims = 1:40)

p <- silhouette_score(seurat_integrated)
pdf(file.path("results", "grid", "k18p12r3_sil.pdf"), height = 12, width = 12)
p
dev.off()

Idents(seurat_integrated) <- "integrated_snn_res.0.4"
seurat_integrated <- RunUMAP(seurat_integrated, dims = 1:40)

p <- silhouette_score(seurat_integrated)
pdf(file.path("results", "grid", "k18p12r4_sil.pdf"), height = 12, width = 12)
p
dev.off()


seurat_integrated <- readRDS(file.path("results", "grid", "seurat_integrated_k18p17.RDS"))
Idents(seurat_integrated) <- "integrated_snn_res.0.3"
seurat_integrated <- RunUMAP(seurat_integrated, dims = 1:40)

p <- silhouette_score(seurat_integrated)
pdf(file.path("results", "grid", "k18p17r3_sil.pdf"), height = 12, width = 12)
p
dev.off()

Idents(seurat_integrated) <- "integrated_snn_res.0.4"
seurat_integrated <- RunUMAP(seurat_integrated, dims = 1:40)

p <- silhouette_score(seurat_integrated)
pdf(file.path("results", "grid", "k18p17r4_sil.pdf"), height = 12, width = 12)
p
dev.off()


seurat_integrated <- readRDS(file.path("results", "grid", "seurat_integrated_k18p17.RDS"))
Idents(seurat_integrated) <- "integrated_snn_res.0.3"
seurat_integrated <- RunUMAP(seurat_integrated, dims = 1:40)

p <- silhouette_score(seurat_integrated)
pdf(file.path("results", "grid", "k18p17r3_sil.pdf"), height = 12, width = 12)
p
dev.off()

Idents(seurat_integrated) <- "integrated_snn_res.0.4"
seurat_integrated <- RunUMAP(seurat_integrated, dims = 1:40)

p <- silhouette_score(seurat_integrated)
pdf(file.path("results", "grid", "k18p17r4_sil.pdf"), height = 12, width = 12)
p
dev.off()

seurat_integrated <- readRDS(file.path("results", "grid", "seurat_integrated_k22.RDS"))
Idents(seurat_integrated) <- "integrated_snn_res.0.3"
seurat_integrated <- RunUMAP(seurat_integrated, dims = 1:40)

p <- silhouette_score(seurat_integrated)
pdf(file.path("results", "grid", "k22r3_sil.pdf"), height = 12, width = 12)
p
dev.off()

Idents(seurat_integrated) <- "integrated_snn_res.0.4"
seurat_integrated <- RunUMAP(seurat_integrated, dims = 1:40)

p <- silhouette_score(seurat_integrated)
pdf(file.path("results", "grid", "k22r4_sil.pdf"), height = 12, width = 12)
p
dev.off()

seurat_integrated <- readRDS(file.path("results", "grid", "seurat_integrated_k22.RDS"))
Idents(seurat_integrated) <- "integrated_snn_res.0.3"
seurat_integrated <- RunUMAP(seurat_integrated, dims = 1:40)

p <- silhouette_score(seurat_integrated)
pdf(file.path("results", "grid", "k22r3_sil.pdf"), height = 12, width = 12)
p
dev.off()

Idents(seurat_integrated) <- "integrated_snn_res.0.4"
seurat_integrated <- RunUMAP(seurat_integrated, dims = 1:40)

p <- silhouette_score(seurat_integrated)
pdf(file.path("results", "grid", "k22r4_sil.pdf"), height = 12, width = 12)
p
dev.off()


seurat_integrated <- readRDS(file.path("results", "grid", "seurat_integrated_k22p12.RDS"))
Idents(seurat_integrated) <- "integrated_snn_res.0.3"
seurat_integrated <- RunUMAP(seurat_integrated, dims = 1:40)

p <- silhouette_score(seurat_integrated)
pdf(file.path("results", "grid", "k22p12r3_sil.pdf"), height = 12, width = 12)
p
dev.off()

Idents(seurat_integrated) <- "integrated_snn_res.0.4"
seurat_integrated <- RunUMAP(seurat_integrated, dims = 1:40)

p <- silhouette_score(seurat_integrated)
pdf(file.path("results", "grid", "k22p12r4_sil.pdf"), height = 12, width = 12)
p
dev.off()

seurat_integrated <- readRDS(file.path("results", "grid", "seurat_integrated_k22p17.RDS"))
Idents(seurat_integrated) <- "integrated_snn_res.0.3"
seurat_integrated <- RunUMAP(seurat_integrated, dims = 1:40)

p <- silhouette_score(seurat_integrated)
pdf(file.path("results", "grid", "k22p17r3_sil.pdf"), height = 12, width = 12)
p
dev.off()

Idents(seurat_integrated) <- "integrated_snn_res.0.4"
seurat_integrated <- RunUMAP(seurat_integrated, dims = 1:40)

p <- silhouette_score(seurat_integrated)
pdf(file.path("results", "grid", "k22p17r4_sil.pdf"), height = 12, width = 12)
p
dev.off()


seurat_integrated <- readRDS(file.path("results", "grid", "seurat_integrated_p12.RDS"))
Idents(seurat_integrated) <- "integrated_snn_res.0.3"
seurat_integrated <- RunUMAP(seurat_integrated, dims = 1:40)

p <- silhouette_score(seurat_integrated)
pdf(file.path("results", "grid", "p12r3_sil.pdf"), height = 12, width = 12)
p
dev.off()

Idents(seurat_integrated) <- "integrated_snn_res.0.4"
seurat_integrated <- RunUMAP(seurat_integrated, dims = 1:40)

p <- silhouette_score(seurat_integrated)
pdf(file.path("results", "grid", "p12r4_sil.pdf"), height = 12, width = 12)
p
dev.off()


