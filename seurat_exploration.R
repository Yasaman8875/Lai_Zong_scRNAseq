library("Seurat")
library("tidyverse")
library("data.table")
library("scProportionTest")
library("future")

options(future.globals.maxSize = 10000 * 1024 ^2)
plan("multiprocess", workers = 12)

seurat_integrated <- readRDS(file.path("results", "r_objects", "seurat_integrated.RDS"))
Idents(seurat_integrated) <- "integrated_snn_res.0.7"
DefaultAssay(seurat_integrated) <- "SCT"
##################
## Marker Genes ##
##################

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
  mutate(avg_fc = (DMSO_avg_logFC + EZH2i_avg_logFC + Combo_avg_logFC + RACi_avg_logFC) / 4) %>% 
  group_by(cluster_id) %>% 
  top_n(n = 10, 
        wt = avg_fc)
write.table(top10, file.path("results", "markers", "top10.csv"), sep=",", col.names=TRUE, 
            row.names=FALSE, quote=FALSE)

## Find all markers
markers <- FindAllMarkers(
  seurat_integrated, assay = "SCT", slot = "data", min.pct = 0.25,
  return.thresh = 0.05, logfc.threshold = log(1.5)
)

saveRDS(markers, file.path("results", "r_objects", "markers.RDS"))
## Load and prepare marker list.

markers <- subset(markers, p_val_adj < 0.05)
ann_markers <- inner_join(x = markers, 
                          y = annotations[, c("gene_name", "description")],
                          by = c("gene" = "gene_name")) %>%
  unique()

# Rearrange the columns to be more intuitive
ann_markers <- ann_markers[ , c(6, 7, 2:4, 1, 5,8)]

# Order the rows by p-adjusted values
ann_markers <- ann_markers %>%
  dplyr::arrange(cluster, p_val_adj)

# Save markers to file
write.csv(ann_markers, 
          file = file.path("results", "markers", "combined_all_markers.csv"), 
          quote = FALSE, 
          row.names = FALSE)

# Extract top 5 markers per cluster
top5 <- ann_markers %>%
  group_by(cluster) %>%
  top_n(n = 5,
        wt = avg_logFC)

write.table(top5, file.path("results", "markers", "top5_combined.csv"), sep=",", 
            col.names=TRUE, row.names=FALSE, quote=FALSE)
## Split the list based on cluster and save results.

ann_markers <- split(ann_markers, ann_markers$cluster)


iwalk(ann_markers, function(x, y) {
  file_name <- file.path("results", "markers", str_c("markers_", y, ".tsv"))
  fwrite(x, file_name, quote = FALSE, sep = "\t", row.names = FALSE, col.names = TRUE)
})

## Cell Cycle Counts
## ----------
## Prepare data for permutation test.

if (!dir.exists(file.path("results", "cell_cycle"))) {
  dir.create(file.path("results", "cell_cycle"))
}

sc_utils_obj <- sc_utils(seurat_integrated)

comparisons <- list(
  c("EZH2i", "DMSO"),
  c("RACi", "DMSO"),
  c("Combo", "DMSO"),
  c("EZH2i", "RACi"),
  c("Combo", "EZH2i"),
  c("Combo", "RACi")
)

## Permutation tests and plotting.

walk(comparisons, function(x) {
  sc_utils_obj <- permutation_test(
    sc_utils_obj, cluster_identity = "Phase",
    sample_1 = x[1], sample_2 = x[2]
  )
  
  p <- permutation_plot(sc_utils_obj, log2FD_threshold = 1)
  
  file_name <- str_c(x[1], "_vs_", x[2], ".pdf")
  pdf(file.path("results", "cell_cycle", file_name), height = 3, width = 8)
  print(p); dev.off()
})

if (!dir.exists(file.path("results", "cluster_counts"))) {
  dir.create(file.path("results", "cluster_counts"))
}

sc_utils_obj <- sc_utils(seurat_integrated)

walk(comparisons, function(x) {
  sc_utils_obj <- permutation_test(
    sc_utils_obj, cluster_identity = "integrated_snn_res.0.7",
    sample_1 = x[1], sample_2 = x[2]
  )
  
  p <- permutation_plot(sc_utils_obj, log2FD_threshold = 1)
  
  file_name <- str_c(x[1], "_vs_", x[2], ".pdf")
  pdf(file.path("results", "cluster_counts", file_name), height = 3, width = 8)
  print(p); dev.off()
})

## Custom Genes Analysis
genes <- c("ALDH1A1", "ALDH1A3", "CD24", "CD44", "PROM1", "BMI1", "SOX2")
# Currently omitting CD117, OCT4, NANOG1
if (!dir.exists(file.path("results", "gene_plots"))) {
  dir.create(file.path("results", "gene_plots"))
}
## Make the feature plot.

DefaultAssay(seurat_integrated) <- "SCT"
p <- FeaturePlot(seurat_integrated, features = genes[1:3], 
                 split.by = "orig.ident", max.cutoff = 3, 
                 cols = c("grey", "red"), pt.size = 0.1)

pdf(file.path("results", "gene_plots", "Feature1.pdf"), height = 16, width = 16)
p
dev.off()

p <- FeaturePlot(seurat_integrated, features = genes[4:7], 
                 split.by = "orig.ident", max.cutoff = 3, 
                 cols = c("grey", "red"), pt.size = 0.1)

pdf(file.path("results", "gene_plots", "Feature2.pdf"), height = 16, width = 16)
p
dev.off()

## Violin Plot
plots <- VlnPlot(seurat_integrated, features = genes, split.by = "orig.ident", 
                 pt.size = 0, combine = FALSE)
p <- wrap_plots(plots = plots, ncol = 3)
pdf(file.path("results", "gene_plots", "Violin.pdf"), height = 16, width = 16)
p
dev.off()