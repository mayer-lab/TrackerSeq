source(file = "/datastore_share/Users/neuhaus/dorsal_ventral_comp/Paper/scripts/lib.R")
library(Seurat)
library(Matrix)
library(ggplot2)
library(dplyr)

out_path1 <- "/data/mayerlab/neuhaus/ferret_proj/final_results/01/"
out_path2 <- "/datastore_share/Users/neuhaus/ferret_proj/final_results/01/"

ferret_seurat <- readRDS(file = "/datastore_share/Users/neuhaus/ferret_proj/data/Ferret_seurat.rds")


########################### SEURAT CLEANING AND CLUSTERING ###########################

## load ref seurat:
ferret_ref <- readRDS(file = "/datastore_share/Users/neuhaus/ferret_proj/data/Mustelaputorius_integrated.rds")

## fix seurat object:
trans_md <- ferret_seurat@meta.data[is.na(ferret_seurat$cloneID), ]
count_assay <- GetAssay(ferret_seurat)
count_mtx <- ferret_seurat@assays$RNA@counts

sum(colnames(count_mtx) == rownames(trans_md))
dim(count_mtx)
## some genes are not unique..
sum(rownames(count_mtx) == make.unique(rownames(count_mtx)))

library(stringr)
str_view_all(rownames(count_mtx), pattern = "[:alnum:]", match = T)
str_view_all(rownames(count_mtx), pattern = "[:space:]", match = T)
str_view_all(rownames(count_mtx), pattern = "[:blank:]", match = T)
str_view_all(rownames(count_mtx), pattern = "[:punct:]", match = T)

new_gene_names <- gsub("_", "-", rownames(count_mtx))
str_view_all(new_gene_names, pattern = "[:punct:]", match = T)
new_gene_names <- make.unique(new_gene_names)
count_mtx@Dimnames[[1]] <- new_gene_names

## only keep genes that also present in ferret ref:
sum(rownames(count_mtx) %in% rownames(ferret_ref)) / nrow(count_mtx)
count_mtx_sub <- count_mtx[rownames(count_mtx) %in% rownames(ferret_ref), ]


ferret_seurat2 <- CreateSeuratObject(
  counts = count_mtx_sub,
  meta.data = trans_md,
  min.cells = 200,
  min.features = 200
)
ferret_seurat2[["RNA"]]@meta.features <- data.frame(row.names = rownames(ferret_seurat2[["RNA"]]))

Idents(ferret_seurat2) <- "orig.ident"

## mitochondrial read fraction:
ferret_seurat2[["percent.mt"]] <- PercentageFeatureSet(ferret_seurat2, pattern = "^MT-")
VlnPlot(ferret_seurat2, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3)
ggsave(paste0(out_path1, "Vln_plot_before_filtering.png"))
ggsave(paste0(out_path2, "Vln_plot_before_filtering.png"))

## filtering:
ferret_seurat2 <- subset(ferret_seurat2, subset = percent.mt < 5 & nFeature_RNA < 4000 & nCount_RNA < 10000 & nFeature_RNA > 1000)
VlnPlot(ferret_seurat2, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3)
ggsave(paste0(out_path1, "Vln_plot_after_filtering.png"))
ggsave(paste0(out_path2, "Vln_plot_after_filtering.png"))

## normalization:
ferret_seurat2 <- NormalizeData(ferret_seurat2, normalization.method = "LogNormalize", scale.factor = 10000)

## variable features:
ferret_seurat2 <- FindVariableFeatures(ferret_seurat2, selection.method = "vst", nfeatures = 2000)
p1 <- VariableFeaturePlot(ferret_seurat2)
p2 <- LabelPoints(plot = p1, points = head(VariableFeatures(ferret_seurat2), 10), repel = TRUE)
p2
ggsave(paste0(out_path1, "Variable_features_plot.png"))
ggsave(paste0(out_path2, "Variable_features_plot.png"))


## regress and scale data:
ferret_seurat2 <- CellCycleScoring(ferret_seurat2, s.features = cc.genes$s.genes, g2m.features = cc.genes$g2m.genes, set.ident = FALSE)
ferret_seurat2$CC.Difference <- ferret_seurat2$S.Score - ferret_seurat2$G2M.Score
## scale data + regress out: nFeature, nCount, percent.mt:
ferret_seurat2 <- ScaleData(ferret_seurat2, vars.to.regress = c("nFeature_RNA", "nCount_RNA", "percent.mt"), verbose = T)


## PCA:
ferret_seurat2 <- RunPCA(ferret_seurat2, features = VariableFeatures(ferret_seurat2), verbose = T)
DimPlot(ferret_seurat2, reduction = "pca", group.by = "orig.ident")
#ggsave(paste0(out_path, "DimPlot_PCA.png"), width = 14, height = 10)

ElbowPlot(ferret_seurat2)
#ggsave(paste0(out_path, "pc_elbow_plot.png"))

## Clustering
ferret_seurat2 <- FindNeighbors(ferret_seurat2, dims = 1:10)
ferret_seurat2 <- FindClusters(ferret_seurat2)

FeaturePlot(ferret_seurat2, features = c("nCount_RNA", "nFeature_RNA"))

ferret_seurat2 <- RunUMAP(ferret_seurat2, dims = 1:10)
DimPlot(ferret_seurat2, reduction = "umap", group.by = "seurat_clusters")
#ggsave(paste0(out_path, "umap_seurat_clusters.png"), width = 14, height = 10)
DimPlot(ferret_seurat2, reduction = "umap", group.by = "orig.ident")
#ggsave(paste0(out_path, "umap_orig_ident.png"), width = 14, height = 10)
DimPlot(ferret_seurat2, reduction = "umap", group.by = "Phase")
#ggsave(paste0(out_path, "umap_phase.png"), width = 14, height = 10)

## save object:
saveRDS(ferret_seurat2, paste0(out_path1, "ferret_seurat2_clean.rds"))


## run doublet detection:
library(DoubletFinder)

## pK Identification (no ground-truth) ---------------------------------------------------------------------------------------
sweep.res.list_ferret <- paramSweep_v3(ferret_seurat2, PCs = 1:10, sct = FALSE)
sweep.stats_ferret <- summarizeSweep(sweep.res.list_ferret, GT = FALSE)
bcmvn_ferret <- find.pK(sweep.stats_ferret)


## Homotypic Doublet Proportion Estimate -------------------------------------------------------------------------------------
annotations <- ferret_seurat2@meta.data$seurat_clusters
homotypic.prop <- modelHomotypic(annotations)
nExp_poi <- round(0.075*nrow(ferret_seurat2@meta.data))  ## Assuming 7.5% doublet formation rate - tailor for your dataset
nExp_poi.adj <- round(nExp_poi*(1-homotypic.prop))

## Run DoubletFinder with varying classification stringencies ----------------------------------------------------------------
ferret_seurat3 <- doubletFinder_v3(ferret_seurat2, PCs = 1:10, pN = 0.25, pK = 0.09, nExp = nExp_poi, reuse.pANN = FALSE, sct = FALSE)
ferret_seurat3 <- doubletFinder_v3(ferret_seurat3, PCs = 1:10, pN = 0.25, pK = 0.09, nExp = nExp_poi.adj, reuse.pANN = "pANN_0.25_0.09_4146", sct = FALSE)

DimPlot(ferret_seurat3, group.by = c("DF.classifications_0.25_0.09_4146","DF.classifications_0.25_0.09_3901"))
table(ferret_seurat3$DF.classifications_0.25_0.09_3901)
table(ferret_seurat3$DF.classifications_0.25_0.09_4146)
DimPlot(ferret_seurat3, group.by = c("DF.classifications_0.25_0.09_3901"))

ferret_seurat_sub <- subset(ferret_seurat3, subset = DF.classifications_0.25_0.09_4146 == "Singlet")
saveRDS(ferret_seurat3, paste0(out_path1, "ferret_seurat3_clean_wDF.rds"))
saveRDS(ferret_seurat_sub, paste0(out_path1, "ferret_seurat_clean_wDF_sub.rds"))


#############################################################################################################################
## re-do scaling and clustering:
VlnPlot(ferret_seurat_sub, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3, group.by = "orig.ident")

ferret_seurat_sub <- FindVariableFeatures(ferret_seurat_sub, selection.method = "vst", nfeatures = 2000)
p1 <- VariableFeaturePlot(ferret_seurat_sub)
p2 <- LabelPoints(plot = p1, points = head(VariableFeatures(ferret_seurat_sub), 10), repel = TRUE)
p2
ggsave(paste0(out_path1, "Variable_features_plot.png"))
ggsave(paste0(out_path2, "Variable_features_plot.png"))

ferret_seurat_sub <- ScaleData(ferret_seurat_sub, vars.to.regress = c("nFeature_RNA", "nCount_RNA", "percent.mt", "CC.Difference"), verbose = T)
## PCA:
ferret_seurat_sub <- RunPCA(ferret_seurat_sub, features = VariableFeatures(ferret_seurat_sub), verbose = T)
DimPlot(ferret_seurat_sub, reduction = "pca", group.by = "orig.ident")
ggsave(paste0(out_path1, "DimPlot_PCA.png"), width = 14, height = 10)
ggsave(paste0(out_path2, "DimPlot_PCA.png"), width = 14, height = 10)

ElbowPlot(ferret_seurat_sub)
ggsave(paste0(out_path1, "pc_elbow_plot.png"))
ggsave(paste0(out_path2, "pc_elbow_plot.png"))

ferret_seurat_sub <- FindNeighbors(ferret_seurat_sub, dims = 1:10)
ferret_seurat_sub <- FindClusters(ferret_seurat_sub)


ferret_seurat_sub <- RunUMAP(ferret_seurat_sub, dims = 1:10)
d1 <- DimPlot(ferret_seurat_sub, reduction = "umap", group.by = "seurat_clusters", label = T) + NoLegend()
d1
ggsave(paste0(out_path1, "umap_seurat_clusters.png"), width = 14, height = 10)
ggsave(paste0(out_path2, "umap_seurat_clusters.png"), width = 14, height = 10)
d2 <- DimPlot(ferret_seurat_sub, reduction = "umap", group.by = "orig.ident")
d2
ggsave(paste0(out_path1, "umap_orig_ident.png"), width = 14, height = 10)
ggsave(paste0(out_path2, "umap_orig_ident.png"), width = 14, height = 10)
d3 <- DimPlot(ferret_seurat_sub, reduction = "umap", group.by = "Phase")
d3
ggsave(paste0(out_path1, "umap_phase.png"), width = 14, height = 10)
ggsave(paste0(out_path2, "umap_phase.png"), width = 14, height = 10)

d1+d2+d3

FeaturePlot(ferret_seurat_sub, features = c("nCount_RNA", "nFeature_RNA"))
ggsave(paste0(out_path1, "count_feature_umap.png"))
ggsave(paste0(out_path2, "count_feature_umap.png"))

marker_gene_vec <- c("HES1", "VIM", "EOMES", "NEUROD6", "S.Score", "G2M.Score")
VlnPlot(ferret_seurat_sub, marker_gene_vec, group.by = "seurat_clusters", ncol = 6)

# broad_vec <- c("EN 1", "RGC G2M", "EN 2", "EN 3", "EN 4", "IPC-EN 1", "RGC S 1", "EN 5", "RGC S 2",
#                "IPC PM", "EN 6", "IPC-EN 2", "EN 7", "EN 8", "IPC S+G2M 1", "IPC S", "IPC S+G2M 2", "RGC-IPC G2M")
# names(broad_vec) <- as.character(seq(0:17))
# ferret_seurat_sub$broad_celltype <- broad_vec[ferret_seurat_sub$seurat_clusters]


FeaturePlot(ferret_seurat_sub, features = marker_gene_vec)

### do finer clustering ###
# ferret_seurat_sub_rc <- FindClusters(ferret_seurat_sub,resolution = 1.3)
# DimPlot(ferret_seurat_sub_rc, reduction = "umap", group.by = c("Phase", "seurat_clusters"))
# 
# VlnPlot(ferret_seurat_sub_rc, marker_gene_vec, group.by = "seurat_clusters", ncol = 3)
# fine_vec <- c(
#   "EN 1", "RGC G2M-S 1", "RGC S 1", "EN 2", "EN 3",
#   "IPC-EN", "IPC-EN 1", "RGC S 2", "EN 4", "EN 5",
#   "EN 6", "EN 7", "IPC-EN 2", "EN 8", "EN 9",
#   "EN 10", "RGC-IPC G2M-S 1", "RGC-IPC G1 1", "EN 11", "RGC G2M-S 2",
#   "RGC-IPC G1 2", "RGC-IPC G2M 1", "RGC G2M 1", "RGC-IPC G2M-S 2", "EN 12"
# )
# names(fine_vec) <- as.character(seq(0:24))
# ferret_seurat_sub_rc$fine_celltype <- fine_vec[ferret_seurat_sub_rc$seurat_clusters]
# 
# DimPlot(ferret_seurat_sub_rc, group.by = "fine_celltype", label = T)
# ggsave(paste0(out_path, "umap_fine_celltype.png"))

saveRDS(ferret_seurat_sub, "/data/mayerlab/neuhaus/ferret_proj/final_results/01/ferret_seurat_clean_wDF_sub.rds")
# saveRDS(ferret_seurat_sub_rc, "/datastore_share/Users/neuhaus/ferret_proj/results/01/ferret_seurat_clean_wDF_wFC_sub.rds")
