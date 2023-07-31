## POST-MITOTIC DIFFERENCES OF RGC OUTPUT ##

library(Seurat)
library(tidyverse)
library(patchwork)
library(UpSetR)
library(dplyr)

## read ferret seurat
ferret_seurat <- readRDS("/data/mayerlab/neuhaus/ferret_proj/final_results/03/ferret_seurat_clean_wLineage.rds")

## add broader clustering:
unique(ferret_seurat$predicted.id)
broad_celltype_vec <- c(
  "EN","RGC-G1","EN","IPC-EN","IPC","RGC-B",
  "IPC-EN","IPC-EN","IPC-EN","RGC-A","IPC","RGC-IPC",
  "RGC-A","RGC-G1","RGC-B","RGC-C","IPC","IPC",
  "RGC-G1"
)

names(broad_celltype_vec) <- unique(ferret_seurat$predicted.id)
broad_celltype_vec
ferret_seurat$broad_trajectory <- broad_celltype_vec[ferret_seurat$predicted.id]


## comapre RGC-A + EN with RGC-B + EN
clone_id_df <- ferret_seurat@meta.data[, c("cloneID", "broad_trajectory")]
clone_id_df %>%
  drop_na(cloneID) %>%
  group_by(cloneID) %>% filter(n() != 1) %>% arrange(broad_trajectory) -> clone_id_df


bct_per_clone <- lapply(unique(clone_id_df$cloneID), function(cid) {
  paste(unique(clone_id_df$broad_trajectory[clone_id_df$cloneID == cid]), collapse = ",")
})
names(bct_per_clone) <- unique(clone_id_df$cloneID)
table(unlist(bct_per_clone))

a_cids <- names(bct_per_clone)[bct_per_clone == "EN,RGC-A"]
b_cids <- names(bct_per_clone)[bct_per_clone == "EN,RGC-B"]

clone_group_vec <- rep(NA, ncol(ferret_seurat))
names(clone_group_vec) <- rownames(ferret_seurat@meta.data)
clone_group_vec[ferret_seurat$cloneID %in% a_cids] <- "EN - RGC A"
clone_group_vec[ferret_seurat$cloneID %in% b_cids] <- "EN - RGC B"

#clone_group_vec[!ferret_seurat$broad_traj_membership %in% c("(1) Upper-layer EN", "General EN")] <- NA
table(clone_group_vec)
ferret_seurat$clone_groups <- clone_group_vec
DimPlot(ferret_seurat, cells.highlight = list("EN-RGC A" = colnames(ferret_seurat)[ferret_seurat$clone_groups == "EN - RGC A"],
                                              "EN-RGC B" = colnames(ferret_seurat)[ferret_seurat$clone_groups == "EN - RGC B"]),
        cols.highlight = c("EN-RGC A" = "blue", "EN-RGC B" = "red"))

out_path <- "/datastore_share/Users/neuhaus/ferret_proj/results/pdf_plots/"
ggsave(paste0(out_path, "DE_clone_groups_umap_wLegend2.png"), width = 9, height = 6)

DimPlot(ferret_seurat, cells.highlight = list("EN-RGC A" = colnames(ferret_seurat)[ferret_seurat$clone_groups == "EN - RGC A"],
                                              "EN-RGC B" = colnames(ferret_seurat)[ferret_seurat$clone_groups == "EN - RGC B"]),
        cols.highlight = c("EN-RGC A" = "blue", "EN-RGC B" = "red")) + NoLegend()
ggsave(paste0(out_path, "DE_clone_groups_umap_woLegend2.png"), width = 9, height = 6)


clone_group_vec[ferret_seurat$broad_trajectory != "EN"] <- NA
ferret_seurat$clone_groups <- clone_group_vec
table(ferret_seurat$clone_groups)

clone_groups_markers <- FindMarkers(ferret_seurat, ident.1 = "EN - RGC A", ident.2 = "EN - RGC B", group.by = "clone_groups", logfc.threshold = 0)
# view results
head(clone_groups_markers)
write.table(clone_groups_markers, file = "/datastore_share/Users/neuhaus/ferret_proj/results/05/EN_RGCA_RGCB_marker_genes.tsv",quote = F)

## plot folcano:
library(ggrepel)
clone_groups_markers$gene <- rownames(clone_groups_markers)
clone_groups_markers$logPval <- -log10(clone_groups_markers$p_val)
p_th <- 1; fc_th <- 0.5
#clone_groups_markers$p_signif <- clone_groups_markers$logPval > p_th
#clone_groups_markers$fc_signif <- abs(clone_groups_markers$avg_log2FC) > fc_th
sign_vec <- sapply(1:nrow(clone_groups_markers), function(idx) {
  if(clone_groups_markers$logPval[idx] > p_th & clone_groups_markers$avg_log2FC[idx] > fc_th) {"EN - RGC A"}
  else if(clone_groups_markers$logPval[idx] > p_th & clone_groups_markers$avg_log2FC[idx] < -fc_th) {"EN - RGC B"}
  else {"not significant"}
})
#clone_groups_markers$signif <- clone_groups_markers$p_signif & clone_groups_markers$fc_signif
clone_groups_markers$significance <- sign_vec
clone_groups_markers$gene[clone_groups_markers$significance == "not significant"] <- NA
ggplot(clone_groups_markers, aes(avg_log2FC, logPval, label = gene, color = significance)) +
  geom_point() +
  geom_text_repel() +
  ylab("-log10(p-value)") +
  xlab("log2FC") +
  theme_bw() +
  scale_color_manual(values = c("red", "blue", "grey"))
ggsave(paste0(out_path, "DE_volcano_plot2.png"), width = 9, height = 6)

## PCA:
EN_seurat_sub <- subset(ferret_seurat, subset = clone_groups %in% c("EN - RGC A", "EN - RGC B"))
dim(EN_seurat_sub)
EN_seurat_sub <- FindVariableFeatures(EN_seurat_sub, selection.method = "vst", nfeatures = 2000)
EN_seurat_sub <- ScaleData(EN_seurat_sub, vars.to.regress = c("nFeature_RNA", "nCount_RNA", "percent.mt", "CC.Difference"), verbose = T)
EN_seurat_sub <- RunPCA(EN_seurat_sub, features = VariableFeatures(EN_seurat_sub), verbose = T, npcs = 10)
DimPlot(EN_seurat_sub, reduction = "pca", group.by = "clone_groups", pt.size = 2, cols = c("red","blue")) +
  ggtitle("PCA embedding of subsetted ENs") + NoLegend()
ggsave(paste0(out_path, "DE_PCA_embedding.png"), width = 9, height = 6)


## look for marker genes gyrus vs. sulcus:
# newborn neurons fig s6:
gyrus_markers <- c("ACTB", "RPS2", "FEZF2", "LHX2", "NFIA", "PCSK1N", "ENSMPUG21393", "NFIB", "NFIX", "PRKX", "TOP2B")
gyrus_markers2<- c("ACTB", "ENSMPUG00000015547", "FEZF2", "LHX2", "NFIA", "PCSK1N", "NFIB", "NFIX", "PRKX", "TOP2B")
sulcus_markers <- c("ENSMPUG19230", "TCF12")
sulcus_markers2<- c("TCF12")

# library(biomaRt)
# ensembl <- useEnsembl(biomart = "ensembl")
# #listDatasets(ensembl)
# mart <- useDataset("mpfuro_gene_ensembl", useMart("ensembl"))
# gene_df <- getBM(attributes = c("ensembl_gene_id", "hgnc_symbol"), mart = mart)
# gene_df[gene_df$hgnc_symbol %in% gyrus_markers,]


clone_groups_markers[rownames(clone_groups_markers) %in% gyrus_markers2, ]
clone_groups_markers[rownames(clone_groups_markers) %in% sulcus_markers2, ]

#paste(clone_groups_markers$gene[clone_groups_markers$signif & clone_groups_markers$avg_log2FC < 0], collapse = ", ")


## PERMUTATION TEST ##
## just shuffle inside labelled EN population:
num_de_genes <- sum(clone_groups_markers$significance != "not significant")
p_th <- 1; fc_th <- 0.5
num_de_genes_sim2 <- c()
sample_vec <- ferret_seurat$clone_groups

for(i in 1:500) {
  sample_vec[!is.na(sample_vec)] <- sample(sample_vec[!is.na(sample_vec)])
  ferret_seurat$clone_groups_sim <- sample_vec
  cgm_df <- FindMarkers(ferret_seurat, ident.1 = "EN - RGC A", ident.2 = "EN - RGC B", group.by = "clone_groups_sim")
  num_de_genes_sim2 <- c(num_de_genes_sim2, sum(-log10(cgm_df$p_val) > p_th & abs(cgm_df$avg_log2FC) > fc_th))
}
hist(num_de_genes_sim2, xlab = "number of DE-genes", main = "Distribution of number of DE-genes", breaks = 20)
abline(v = num_de_genes, col = "red")

saveRDS(num_de_genes_sim2, file = "/datastore_share/Users/neuhaus/ferret_proj/results/05/num_de_genes_simulated500.rds")
#num_de_genes_sim2 <- readRDS(file = "/datastore_share/Users/neuhaus/ferret_proj/results/05/num_de_genes_simulated500.rds")

num_de_genes_df <- data.frame(
  "num_de_genes" = num_de_genes_sim2
)
ggplot(num_de_genes_df, aes(num_de_genes)) +
  geom_histogram() + 
  theme_bw() +
  xlab("number of DEGs") + ggtitle("Distribution of the number of DEGs") +
  geom_vline(xintercept = num_de_genes, color = "red")
ggsave(paste0(out_path, "DE_num_de_genes_sim2.png"), width = 9, height = 6)
