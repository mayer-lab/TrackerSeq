library(Seurat)
library(ggplot2)
library(dplyr)

ferret_seurat_sub <- readRDS("/data/mayerlab/neuhaus/ferret_proj/final_results/01/ferret_seurat_clean_wDF_sub.rds")
ferret_ref <- readRDS(file = "/datastore_share/Users/neuhaus/ferret_proj/data/Mustelaputorius_integrated.rds")

ferret_ref <- FindVariableFeatures(ferret_ref, selection.method = "vst", nfeatures = 2000)

## exclude interneuron clusters from ref for easier label transfer:
ferret_ref_sub <- subset(ferret_ref, subset = cell_identity %in% 
                           unique(ferret_ref$cell_identity)[!unique(ferret_ref$cell_identity) %in% 
                                                              c("(6) Immature MGE IN SST-","(15) Immature CGE IN","(20) Immature MGE IN SST+NPY/RELN", "(25) OPC", "(26) EC", "(28) Microglia")]
                        )

out_path1 <- "/data/mayerlab/neuhaus/ferret_proj/final_results/02/"
out_path2 <- "/datastore_share/Users/neuhaus/ferret_proj/final_results/02/"

########################### CELL TYPE CLASSIFICATION WITH LABEL TRANSFER ###########################

DimPlot(ferret_ref_sub, reduction = "umap", group.by = "cell_identity", label = T) + NoLegend()

ferret_anchors <- FindTransferAnchors(reference = ferret_ref_sub, query = ferret_seurat_sub)
# Performing PCA on the provided reference using 1597 features as input.
# Projecting cell embeddings
# Finding neighborhoods
# Finding anchors
# Found 15328 anchors
# Filtering anchors
# Retained 7087 anchors

predictions <- TransferData(anchorset = ferret_anchors, refdata = ferret_ref_sub$cell_identity)

ferret_seurat_sub <- AddMetaData(ferret_seurat_sub, metadata = predictions)

table(ferret_seurat_sub$predicted.id)
# (0) IPC-Newborn EN         (1) Upper-layer EN                (10) IPC G1            (11) RGC S/G2/M            (13) IPC S/G2/M 
# 7710                      17141                       1521                       1985                       1506 
# (14) IPC S/G2/M                (16) RGC G1            (17) RGC S/G2/M               (18) RGC-IPC (19) IPC S/G2/M-Newborn EN 
# 522                        180                       1394                       1205                       1004 
# (2) IPC-Newborn EN                (21) Neuron            (22) RGC S/G2/M            (23) RGC S/G2/M            (24) IPC S/G2/M 
# 2897                         81                        667                        513                        264 
# (3) IPC-Newborn EN                 (4) RGC G1             (7) General EN                 (8) RGC G1             (9) RGC S/G2/M 
# 2296                       2757                       4564                       1272                       1653


DimPlot(ferret_seurat_sub, reduction = "umap", group.by = "predicted.id", label = F)
ggsave(paste0(out_path1, "umap_predicted_id_woLabel2.png"), width = 14, height = 10)
ggsave(paste0(out_path2, "umap_predicted_id_woLabel2.png"), width = 14, height = 10)

DimPlot(ferret_seurat_sub, reduction = "umap", group.by = "predicted.id", label = T)
ggsave(paste0(out_path1, "umap_predicted_id_wLabel2.png"), width = 14, height = 10)
ggsave(paste0(out_path2, "umap_predicted_id_wLabel2.png"), width = 14, height = 10)


saveRDS(ferret_seurat_sub, "/data/mayerlab/neuhaus/ferret_proj/final_results/02/ferret_seurat_clean_wLabel.rds")

## try to reconstruct cluster identity with marker genes from manuscript ##
# traj 1: 
# 23: CDC20,CCNB2,CCNB1,ASPM,Pttg1
# 11: BIRC5, Pttg1, CDCA3, HMGN3, VIM

# traj 2:
# 22: CENPF, NUSAP1, TOP2A, ASPM, CCNB1
# 9: ENSMPUG00000018619, HIST3H2A, 2810417H13Rik, TOP2A, DUT

# RGC G1: HES1, FOS, FOSB, VIM, PON2

# RGC-IPC: RGS16, ILKAP, HES1, RPS27L, VIM, CCND1

# IPC G1: NEUROD4, ELAVL4, EOMES, ILKAP

# IPC S/G2/M: DLGAP5, HMMR, TOP2A, CLSPN, NEUROD4, PCNA

# (19) IPC S/G2/M-Newborn EN: BIRC5, Pttg1, NEUROD4, SPC24, CDCA3

# IPC-Newborn EN: EOMES, NRN1, NEUROD1, NEUROD2, CLIC1, DCX, ELAVL4

# EN: NEUROD6, CSRP2, GAP43, NTM, Tubb3

# 23,11:
FeaturePlot(ferret_seurat, features = c("CDC20","CCNB2","ASPM","Pttg1",
                                        "BIRC5", "Pttg1", "CDCA3", "VIM"), ncol = 4)

# 22,9:
FeaturePlot(ferret_seurat, features = c("CENPF", "NUSAP1", "TOP2A", "ASPM", "CCNB1",
                                        "ENSMPUG00000018619", "HIST3H2A", "2810417H13Rik", "TOP2A", "DUT"), ncol = 5)

#RGC G1:
FeaturePlot(ferret_seurat, features = c( "HES1", "FOS", "FOSB", "VIM", "PON2"))

# RGC-IPC: 
FeaturePlot(ferret_seurat, features = c("RGS16", "ILKAP", "HES1", "RPS27L", "VIM", "CCND1"))

# IPC G1:
FeaturePlot(ferret_seurat, features = c( "NEUROD4", "ELAVL4", "EOMES", "ILKAP"))

#  IPC S/G2/M: 
FeaturePlot(ferret_seurat, features = c("DLGAP5", "HMMR", "TOP2A", "CLSPN", "NEUROD4", "PCNA"))

# (19) IPC S/G2/M-Newborn EN: 
FeaturePlot(ferret_seurat, features = c("BIRC5", "Pttg1", "NEUROD4", "SPC24", "CDCA3"))

# IPC-Newborn EN: 
FeaturePlot(ferret_seurat, features = c("EOMES", "NRN1", "NEUROD1", "NEUROD2", "CLIC1", "DCX", "ELAVL4"))

# EN: 
FeaturePlot(ferret_seurat, features = c("NEUROD6", "CSRP2", "GAP43", "NTM", "Tubb3"))


## fig s11:
traj_genes <- c("PTTG1","TOX3","LGALS1","SLC1A3","ZFP36L1","SPARC","TOP2A","NSL1","SPC25")
FeaturePlot(ferret_seurat_sub, features = traj_genes, ncol = 3)



