library("Seurat")
library("tidyverse")
library("data.table")
library("scProportionTest")
library("future")
library("wesanderson")
library(patchwork)

options(future.globals.maxSize = 10000 * 1024 ^2)
plan("multiprocess", workers = 12)

seurat_integrated <- readRDS(file.path("results", "r_objects", "seurat_integrated.RDS"))
Idents(seurat_integrated) <- "integrated_snn_res.0.3"
DefaultAssay(seurat_integrated) <- "SCT"
##################
## Marker Genes ##
##################

if (!dir.exists(file.path("results", "markers"))) {
  dir.create(file.path("results", "markers"))
}
annotations <- read.csv("annotation.csv")

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


# Prints the cluster counts separated per sample.
table(seurat_integrated@meta.data$integrated_snn_res.0.3, seurat_integrated@meta.data$orig.ident)

## Differential Expression Analysis
seurat_integrated$celltype.stim <- paste(Idents(seurat_integrated), seurat_integrated$orig.ident, sep = "_")
seurat_integrated$celltype <- Idents(seurat_integrated)
Idents(seurat_integrated) <- "celltype.stim"

for (i in 1:11) {
  write.table(FindMarkers(seurat_integrated, ident.1 = paste(i,"_EZH2i",sep = ""), 
                          ident.2 = paste(i,"_DMSO",sep = ""), verbose = FALSE, 
                          assay = "SCT", slot = "data"), file.path("results", "markers", 
                          paste(i, "_EZH2i_DMSO.tsv", sep="")), sep = "\t",  
                          col.names=TRUE, row.names=TRUE)
  write.table(FindMarkers(seurat_integrated, ident.1 = paste(i,"_RACi",sep = ""), 
                          ident.2 = paste(i,"_DMSO",sep = ""), verbose = FALSE, 
                          assay = "SCT", slot = "data"), file.path("results", "markers", 
                          paste(i, "_RACi_DMSO.tsv", sep="")), sep = "\t",  
                          col.names=TRUE, row.names=TRUE)
  write.table(FindMarkers(seurat_integrated, ident.1 = paste(i,"_Combo",sep = ""), 
                          ident.2 = paste(i,"_DMSO",sep = ""), verbose = FALSE, 
                          assay = "SCT", slot = "data"), file.path("results", "markers", 
                          paste(i, "_Combo_DMSO.tsv", sep="")), sep = "\t",  
                          col.names=TRUE, row.names=TRUE)
  write.table(FindMarkers(seurat_integrated, ident.1 = paste(i,"_Combo",sep = ""),
			  ident.2 = paste(i,"_RACi",sep = ""), verbose = FALSE,
			  assay = "SCT", slot = "data"), file.path("results", "markers",
			  paste(i, "_Combo_RACi.tsv", sep="")), sep = "\t",
	                  col.names=TRUE, row.names=TRUE)
}
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
p <- wrap_plots(plots = plots, ncol = 2)
pdf(file.path("results", "gene_plots", "Violin.pdf"), height = 16, width = 16)
p
dev.off()

## Plots for cluster 2
genes <- c("TM4SF1", "NMU", "UPK1B")
p <- FeaturePlot(seurat_integrated, features = genes, 
                 split.by = "orig.ident", max.cutoff = 3, 
                 cols = c("grey", "red"), pt.size = 0.1)

pdf(file.path("results", "gene_plots", "Feature_2.pdf"), height = 16, width = 16)
p
dev.off()

plots <- VlnPlot(seurat_integrated, features = genes, split.by = "orig.ident", 
                 pt.size = 0, combine = FALSE)
p <- wrap_plots(plots = plots, ncol = 1)
pdf(file.path("results", "gene_plots", "Violin_2.pdf"), height = 16, width = 16)
p
dev.off()

## Plots for cluster 7
genes <- c("CXCL1", "CXCL2", "CXCL3", "CXCL8")
p <- FeaturePlot(seurat_integrated, features = genes, 
                 split.by = "orig.ident", max.cutoff = 3, 
                 cols = c("grey", "red"), pt.size = 0.1)

pdf(file.path("results", "gene_plots", "Feature_7.pdf"), height = 16, width = 16)
p
dev.off()

plots <- VlnPlot(seurat_integrated, features = genes, split.by = "orig.ident", 
                 pt.size = 0, combine = FALSE)
p <- wrap_plots(plots = plots, ncol = 2)
pdf(file.path("results", "gene_plots", "Violin_7.pdf"), height = 12, width = 20)
p
dev.off()

## Plots for cluster 8
genes <- c("FLG", "CDH6", "TMEM100", "CD24", "CTH")
p <- FeaturePlot(seurat_integrated, features = genes, 
                 split.by = "orig.ident", max.cutoff = 3, 
                 cols = c("grey", "red"), pt.size = 0.1)

pdf(file.path("results", "gene_plots", "Feature_8.pdf"), height = 16, width = 16)
p
dev.off()

plots <- VlnPlot(seurat_integrated, features = genes, split.by = "orig.ident", 
                 pt.size = 0, combine = FALSE)
p <- wrap_plots(plots = plots, ncol = 2)
pdf(file.path("results", "gene_plots", "Violin_8.pdf"), height = 12, width = 20)
p
dev.off()

## Misc 
genes <- c("CDK1", "CCNB1", "HIST2H2AC", "UCA1", "PLCG2")
p <- FeaturePlot(seurat_integrated, features = genes, 
                 split.by = "orig.ident", max.cutoff = 3, 
                 cols = c("grey", "red"), pt.size = 0.1)

pdf(file.path("results", "gene_plots", "FeatureMisc.pdf"), height = 16, width = 16)
p
dev.off()

plots <- VlnPlot(seurat_integrated, features = genes, split.by = "orig.ident", 
                 pt.size = 0, combine = FALSE)
p <- wrap_plots(plots = plots, ncol = 2)
pdf(file.path("results", "gene_plots", "ViolinMisc.pdf"), height = 12, width = 20)
p
dev.off()

# Epithelial analysis
# Early Secretory Epithelial 
genes <- c("PAX8", "THY1")
p <- FeaturePlot(seurat_integrated, features = genes, 
                 split.by = "orig.ident", max.cutoff = 3, 
                 cols = c("grey", "red"), pt.size = 0.1)

pdf(file.path("results", "gene_plots", "Feature_EarlySecEp.pdf"), height = 12, width = 20)
p
dev.off()

## Violin Plot
plots <- VlnPlot(seurat_integrated, features = genes, split.by = "orig.ident", 
                 pt.size = 0, combine = FALSE)
p <- wrap_plots(plots = plots, ncol = 1)
pdf(file.path("results", "gene_plots", "Violin_EarlySecEp.pdf"), height = 20, width = 12)
p
dev.off()

p <- FeaturePlot(seurat_integrated, features = "OVGP1", 
                 split.by = "orig.ident", max.cutoff = 3, 
                 cols = c("grey", "red"), pt.size = 0.1)

pdf(file.path("results", "gene_plots", "Feature_SecEp.pdf"), height = 12, width = 20)
p
dev.off()

## Violin Plot
plots <- VlnPlot(seurat_integrated, features = "OVGP1", split.by = "orig.ident", 
                 pt.size = 0, combine = FALSE)
p <- wrap_plots(plots = plots, ncol = 1)
pdf(file.path("results", "gene_plots", "Violin_SecEp.pdf"), height = 12, width = 12)
p
dev.off()
# Ciliated
genes <- c("FOXJ1", "RUNX3")
p <- FeaturePlot(seurat_integrated, features = genes, 
                 split.by = "orig.ident", max.cutoff = 3, 
                 cols = c("grey", "red"), pt.size = 0.1)

pdf(file.path("results", "gene_plots", "Feature_CilEp.pdf"), height = 10, width = 20)
p
dev.off()

## Violin Plot
plots <- VlnPlot(seurat_integrated, features = genes, split.by = "orig.ident", 
                 pt.size = 0, combine = FALSE)
p <- wrap_plots(plots = plots, ncol = 1)
pdf(file.path("results", "gene_plots", "Violin_CilEp.pdf"), height = 12, width = 12)
p
dev.off()

# Keratin

genes <- c("KRT7", "KRT8", "KRT18")
p <- FeaturePlot(seurat_integrated, features = genes, 
                 split.by = "orig.ident", max.cutoff = 3, 
                 cols = c("grey", "red"), pt.size = 0.1)

pdf(file.path("results", "gene_plots", "Feature_KrtEp.pdf"), height = 18, width = 20)
p
dev.off()

## Violin Plot
plots <- VlnPlot(seurat_integrated, features = genes, split.by = "orig.ident", 
                 pt.size = 0, combine = FALSE)
p <- wrap_plots(plots = plots, ncol = 1)
pdf(file.path("results", "gene_plots", "Violin_KrtEp.pdf"), height = 24, width = 12)
p
dev.off()

# Misc/Unclassified

genes <- c("EPCAM", "KRT10")
p <- FeaturePlot(seurat_integrated, features = genes, 
                 split.by = "orig.ident", max.cutoff = 3, 
                 cols = c("grey", "red"), pt.size = 0.1)

pdf(file.path("results", "gene_plots", "Feature_MiscEp.pdf"), height = 15, width = 20)
p
dev.off()

## Violin Plot
plots <- VlnPlot(seurat_integrated, features = genes, split.by = "orig.ident", 
                 pt.size = 0, combine = FALSE)
p <- wrap_plots(plots = plots, ncol = 1)
pdf(file.path("results", "gene_plots", "Violin_MiscEp.pdf"), height = 16, width = 12)
p
dev.off()

# Label clusters
seurat_integrated <- RenameIdents(object = seurat_integrated, 
                                  "1" = "Misc. 1",
                                  "2" = "CDD",
                                  "3" = "Misc. 2",
                                  "4" = "Cell Cycle",
                                  "5" = "Inflammatory Response",
                                  "6" = "Misc. 3",
                                  "7" = "Misc. 4",
                                  "8" = "Cell Division",
                                  "9" = "GTSE1 enriched",
                                  "10" = "Misc. 5",
                                  "11" = "Doublets")
# Plot the UMAP
p <- DimPlot(object = seurat_integrated, 
            group.by = "ident", 
            split.by = "orig.ident", 
            ncol = 2,
            reduction = "umap")

pdf(file.path("results", "clustering", "clusters.pdf"), height = 12, width = 10)
p
dev.off()

# Integrated clusters
p <- DimPlot(seurat_integrated, 
             reduction = "umap",
             label = FALSE,
             label.size = 3,
             repel = TRUE)
pdf(file.path("results", "clustering", "integrated_clusters.pdf"), height = 10, width = 10)
p
dev.off()

## Cell Cycle Counts
## ----------
## Prepare data for permutation test.

if (!dir.exists(file.path("results", "cell_cycle"))) {
    dir.create(file.path("results", "cell_cycle"))
}

sc_utils_obj <- sc_utils(seurat_integrated)

comparisons <- list(
    c("DMSO", "GSK126"),
    c("DMSO", "NSC23766"),
    c("DMSO", "Combo"),
    c("GSK126", "NSC23766"),
    c("Combo", "GSK126"),
    c("Combo", "NSC23766")
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
seurat_integrated$custom_clusters <- Idents(seurat_integrated)
Idents(seurat_integrated) <- "custom_clusters"

sc_utils_obj <- sc_utils(seurat_integrated)

walk(comparisons, function(x) {
    sc_utils_obj <- permutation_test(
        sc_utils_obj, cluster_identity = "custom_clusters",
        sample_1 = x[1], sample_2 = x[2]
    )
    
    p <- permutation_plot(sc_utils_obj, log2FD_threshold = 1)
    
    file_name <- str_c(x[1], "_vs_", x[2], ".pdf")
    pdf(file.path("results", "cluster_counts", file_name), height = 3, width = 8)
    print(p); dev.off()
})

# Final feature/violin plots for publication

p <- FeaturePlot(seurat_integrated, features = "ALDH1A1", 
                 split.by = "orig.ident", max.cutoff = 3, 
                 cols = c("grey", "red"), pt.size = 0.1)

pdf(file.path("results", "gene_plots", "Feature_ALDH1A1.pdf"), height = 6, width = 15)
p
dev.off()

p <- FeaturePlot(seurat_integrated, features = "SOX2", 
                 split.by = "orig.ident", max.cutoff = 3, 
                 cols = c("grey", "red"), pt.size = 0.1)

pdf(file.path("results", "gene_plots", "Feature_SOX2.pdf"), height = 6, width = 15)
p
dev.off()

p <- FeaturePlot(seurat_integrated, features = "CD24", 
                 split.by = "orig.ident", max.cutoff = 3, 
                 cols = c("grey", "red"), pt.size = 0.1)

pdf(file.path("results", "gene_plots", "Feature_CD24.pdf"), height = 6, width = 15)
p
dev.off()

genes <- c("PAX8", "OVGP1")

p <- FeaturePlot(seurat_integrated, features = genes, 
                 split.by = "orig.ident", max.cutoff = 3, 
                 cols = c("grey", "red"), pt.size = 0.1)

pdf(file.path("results", "gene_plots", "Feature_PAX8_OVGP1.pdf"), height = 12, width = 15)
p
dev.off()

p <- FeaturePlot(seurat_integrated, features = genes, 
                 split.by = "orig.ident", max.cutoff = 3, 
                 cols = c("grey", "red"), pt.size = 0.1)

pdf(file.path("results", "gene_plots", "Feature_OVGP1_OVGP1.pdf"), height = 12, width = 15)
p
dev.off()

# Violin Plots

plots <- VlnPlot(seurat_integrated, features = "ALDH1A1", split.by = "orig.ident", 
                 pt.size = 0, combine = FALSE)
p <- wrap_plots(plots = plots, ncol = 1)
pdf(file.path("results", "gene_plots", "Violin_ALDH1A1.pdf"), height = 8, width = 12)
p
dev.off()

plots <- VlnPlot(seurat_integrated, features = "SOX2", split.by = "orig.ident", 
                 pt.size = 0, combine = FALSE)
p <- wrap_plots(plots = plots, ncol = 1)
pdf(file.path("results", "gene_plots", "Violin_SOX2.pdf"), height = 8, width = 12)
p
dev.off()

plots <- VlnPlot(seurat_integrated, features = "CD24", split.by = "orig.ident", 
                 pt.size = 0, combine = FALSE)
p <- wrap_plots(plots = plots, ncol = 1)
pdf(file.path("results", "gene_plots", "Violin_CD24.pdf"), height = 8, width = 12)
p
dev.off()

plots <- VlnPlot(seurat_integrated, features = "PAX8", split.by = "orig.ident", 
                 pt.size = 0, combine = FALSE)
p <- wrap_plots(plots = plots, ncol = 1)
pdf(file.path("results", "gene_plots", "Violin_PAX8.pdf"), height = 8, width = 12)
p
dev.off()

plots <- VlnPlot(seurat_integrated, features = "OVGP1", split.by = "orig.ident", 
                 pt.size = 0, combine = FALSE)
p <- wrap_plots(plots = plots, ncol = 1)
pdf(file.path("results", "gene_plots", "Violin_OVGP1.pdf"), height = 8, width = 12)
p
dev.off()

p <- FeaturePlot(seurat_integrated, features = "SIRT7", 
                 split.by = "orig.ident", max.cutoff = 3, 
                 cols = c("grey", "red"), pt.size = 0.1)

pdf(file.path("results", "gene_plots", "Feature_SIRT7.pdf"), height = 12, width = 20)
p
dev.off()

plots <- VlnPlot(seurat_integrated, features = "SIRT7", split.by = "orig.ident", 
                 pt.size = 0, combine = FALSE)
p <- wrap_plots(plots = plots, ncol = 1)
pdf(file.path("results", "gene_plots", "Violin_SIRT7.pdf"), height = 12, width = 12)
p
dev.off()

# Consider SIRT7 in cluster 2 specifically
seurat_integrated <- subset(seurat_integrated, subset = integrated_snn_res.0.3 == "2")
p <- FeaturePlot(seurat_integrated, features = "SIRT7", 
                 split.by = "orig.ident", max.cutoff = 3, 
                 cols = c("grey", "red"), pt.size = 0.1)

pdf(file.path("results", "gene_plots", "Cluster2_SIRT7.pdf"), height = 12, width = 20)
p
dev.off()


## Cell Cycle

if (!dir.exists(file.path("results", "cell_cycle"))) {
	dir.create(file.path("results", "cell_cycle"))
}

cell_cycle_palette <- wes_palette("Zissou1", 3, type = "continuous") %>%
	as.character

## Dim plots of cell cycle phase.

p <- DimPlot(
	seurat_integrated, group.by = "Phase", split.by = "orig.ident",
	ncol = 2, cols = cell_cycle_palette
)

pdf(file.path("results", "cell_cycle", "cell_cycle_dimplot.pdf"), height = 10, width = 10)
p
dev.off()

