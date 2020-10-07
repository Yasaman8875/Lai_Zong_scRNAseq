library(Seurat)
library("tidyverse")
library("clustree")
library("future")
library("unixtools")
library(ggplot2)
library(patchwork)
library("scProportionTest")

options(future.globals.maxSize = 10000 * 1024 ^2)
plan("multiprocess", workers = 4)

if (!dir.exists("results")) dir.create("results")

samples_10x <- list(DMSO = file.path("DMSO", "outs", "filtered_feature_bc_matrix"),
		    EZH2 = file.path("EZH2", "outs", "filtered_feature_bc_matrix"),
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
mt_data <- map(seurat_obj, function(x) {
	meta_data <- x@meta.data
	meta_data <- as_tibble(meta_data, .name_repair = "unique")
        meta_data$cells <- rownames(meta_data) 
        meta_data <- meta_data %>% rename(nUmi = nCount_RNA, nGene = nFeature_RNA)
	return(meta_data)
})

mt_data$DMSO$sample <- "DMSO"
mt_data$EZH2$sample <- "EZH2"
mt_data$RACi$sample <- "RACi"
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

