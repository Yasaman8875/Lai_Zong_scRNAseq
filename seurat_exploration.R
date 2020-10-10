library("Seurat")
library("tidyverse")
library("data.table")
library("scProportionTest")
library("future")

options(future.globals.maxSize = 10000 * 1024 ^2)
plan("multiprocess", workers = 12)

seurat_integrated <- readRDS(file.path("results", "r_objects", "seurat_integrated.RDS"))

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
  mutate(avg_fc = (DMSO_avg_logFC + EZH2i_avg_logFC + Combo_avg_logFC + RACi_avg_logFC) / 4) %>% 
  group_by(cluster_id) %>% 
  top_n(n = 10, 
        wt = avg_fc)
write.table(top10, file.path("results", "markers", "top10.csv"), sep=",", col.names=TRUE, 
            row.names=FALSE, quote=FALSE)

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
