library(Seurat)
library(tidyverse)
library(patchwork)
library(UpSetR)
library(dplyr)

## read ferret seurat with label transfer:
ferret_seurat <- readRDS("/data/mayerlab/neuhaus/ferret_proj/final_results/02/ferret_seurat_clean_wLabel.rds")
#ferret_seurat <- readRDS("/datastore_share/Users/neuhaus/ferret_proj/results/03/ferret_seurat_clean_wLabel_wLineage.rds")

out_path1 <- "/data/mayerlab/neuhaus/ferret_proj/final_results/03/"
out_path2 <- "/datastore_share/Users/neuhaus/ferret_proj/final_results/03/"

########################### LINEAGE ANALYSIS ###########################

clones <- read.table("/data/mayerlab/neuhaus/ferret_proj/final_results/00/ED221124_cloneIDs_nreads10_numi9.csv", h=T, row.names = 1, sep = ",")
table(table(clones$cloneID_filtered))

## add clone-id to seurat:
orig_cell_ids <- sapply(colnames(ferret_seurat), function(el) {strsplit(el, "_")[[1]][2]})

## create clone ID vector:
create_cloneID_vec <- function(clones_sub, orig_cell_ids, ferret_seurat, filtered = FALSE) {
   clone_id_vec <- rep(NA, length(orig_cell_ids))
   names(clone_id_vec) <- orig_cell_ids
   if(filtered) {
      clone_id_vec[clones_sub$cellbc] <- clones_sub$cloneID_filtered
   } else {
      clone_id_vec[clones_sub$cellbc] <- clones_sub$cloneID
   }
   names(clone_id_vec) <- colnames(ferret_seurat)
   return(clone_id_vec)
}
ferret_seurat$cloneID <- create_cloneID_vec(clones, orig_cell_ids, ferret_seurat, filtered = TRUE)


## presentation stuff:
DimPlot(ferret_seurat, cells.highlight = colnames(ferret_seurat)[!is.na(ferret_seurat$cloneID)])
ggsave(paste0(out_path1, "cells_with_cloneIDs_dimplot.png"))
ggsave(paste0(out_path2, "cells_with_cloneIDs_dimplot.png"))


cells_per_cluster <- data.frame(
   "predicted_id" = ferret_seurat$predicted.id,
   "clone_id" = ferret_seurat$cloneID
)
cells_per_cluster <- cells_per_cluster[!is.na(cells_per_cluster$clone_id), ]
cells_per_pid <- sapply(unique(cells_per_cluster$predicted_id), function(pid) {
   sum(cells_per_cluster$predicted_id == pid)
})
names(cells_per_pid) <- unique(cells_per_cluster$predicted_id)

cells_per_pid_df <- data.frame("predicted_id" = names(cells_per_pid), "count" = cells_per_pid)
ggplot(cells_per_pid_df, aes(predicted_id, count)) +
   geom_bar(stat = "identity") +
   theme_bw() + theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))
ggsave(paste0(out_path1, "cells_per_predictedID_barplot.png"))
ggsave(paste0(out_path2, "cells_per_predictedID_barplot.png"))


saveRDS(ferret_seurat, "/data/mayerlab/neuhaus/ferret_proj/final_results/03/ferret_seurat_clean_wLineage.rds")


############################# UPSET PLOT ############################# 

pool_upset <- FetchData(ferret_seurat, c("cloneID", "predicted.id"))
pool_upset <- pool_upset %>% drop_na(cloneID)
pool_upset <- pool_upset %>% ungroup() %>% group_by(cloneID) %>% filter(n() != 1)

png(paste0(out_path2, "upset_plot_all_predictedIDs.png"), width = 18, height = 10)
upset(data = fromList(split(pool_upset$cloneID, pool_upset$predicted.id)),
      nintersects = 148,
      nsets = 100,
      order.by = "freq",
      decreasing = T,
      mb.ratio = c(0.3,0.5)
)
dev.off()

## upset for traj1:
traj1_clusters <- c("(7) General EN", "(8) RGC G1","(1) Upper-layer EN","(2) IPC-Newborn EN", "(19) IPC S/G2/M-Newborn EN", "(0) IPC-Newborn EN","(3) IPC-Newborn EN" ,"(11) RGC S/G2/M","(18) RGC-IPC","(23) RGC S/G2/M", "(17) RGC S/G2/M","(10) IPC G1")
traj2_clusters <- c("(7) General EN", "(8) RGC G1","(1) Upper-layer EN","(2) IPC-Newborn EN", "(19) IPC S/G2/M-Newborn EN", "(0) IPC-Newborn EN","(3) IPC-Newborn EN" ,"(9) RGC S/G2/M","(18) RGC-IPC","(22) RGC S/G2/M", "(17) RGC S/G2/M","(10) IPC G1", "(4) RGC G1", "(16) RGC G1")
traj3_clusters <- c("(22) RGC S/G2/M", "(9) RGC S/G2/M", "(13) IPC S/G2/M", "(24) IPC S/G2/M", "(14) IPC S/G2/M", "(7) General EN", "(10) IPC G1", "(19) IPC S/G2/M-Newborn EN", "(2) IPC-Newborn EN", "(1) Upper-layer EN", "(0) IPC-Newborn EN","(3) IPC-Newborn EN")

pool_upset2_t1 <- pool_upset[pool_upset$predicted.id %in% traj1_clusters, ]
pool_upset2_t2 <- pool_upset[pool_upset$predicted.id %in% traj2_clusters, ]
pool_upset2_t3 <- pool_upset[pool_upset$predicted.id %in% traj3_clusters, ]

png(paste0(out_path2, "upset_plot_trajectory1_predictedID.png"), width = 18, height = 10)
upset(data = fromList(split(pool_upset2_t1$cloneID, pool_upset2_t1$predicted.id)),
      nintersects = 100,
      nsets = 100,
      order.by = "freq",
      decreasing = T,
      mb.ratio = c(0.3,0.5)
)
dev.off()

png(paste0(out_path2, "upset_plot_trajectory2_predictedID.png"), width = 18, height = 10)
upset(data = fromList(split(pool_upset2_t2$cloneID, pool_upset2_t2$predicted.id)),
      nintersects = 100,
      nsets = 100,
      order.by = "freq",
      decreasing = T,
      mb.ratio = c(0.3,0.5)
)
dev.off()

png(paste0(out_path2, "upset_plot_trajectory3_predictedID.png"), width = 18, height = 10)
upset(data = fromList(split(pool_upset2_t3$cloneID, pool_upset2_t3$predicted.id)),
      nintersects = 100,
      nsets = 100,
      order.by = "freq",
      decreasing = T,
      mb.ratio = c(0.3,0.5)
)
dev.off()

###############################################################################################################
#########################        AFTER LINEAGE COUPLING                  ######################################
###############################################################################################################

library(data.table)
export_for_cytoscape <- function(cor_mtx) {
   lc_tl_broad <- as.matrix(cor_mtx)
   lc_tl_broad[lower.tri(lc_tl_broad, diag = T)] <- NA
   lc_tl_broad <- as.data.table(lc_tl_broad)
   lc_tl_broad <- lc_tl_broad[, "source" := rownames(cor_mtx)]
   lc_tl_broad %>%
      data.table::melt(value.name = "correlation", id.var = "source", variable.name = "target") -> lc_tl_long
   lc_tl_long <- lc_tl_long[!is.na(correlation)]
   return(lc_tl_long)
}



lc_tl_cor_mtx <- read.table("/data/mayerlab/neuhaus/ferret_proj/final_results/04/lineage_coupling_scores_correlation_matrix.csv", h = T, row.names = 1, sep = ",")
colnames(lc_tl_cor_mtx) <- rownames(lc_tl_cor_mtx)
library(pheatmap)
png(paste0(out_path2, "lineage_correlation_heatmap.png"), width = 1100, height = 1100)
pheatmap(lc_tl_cor_mtx, cluster_rows = T, cluster_cols = T, display_numbers = T, breaks = seq(-1,1,length.out = 100))
dev.off()

lc_tl_cor_long <- export_for_cytoscape(lc_tl_cor_mtx)
write.table(lc_tl_cor_long, paste0(out_path2, "lc_tl_cor_matrix_cytoscape.csv"), sep = ",", quote = F, row.names = F)


lc_tl_z_mtx <- read.table("/data/mayerlab/neuhaus/ferret_proj/final_results/04/lineage_coupling_scores_matrix.csv", h = T, row.names = 1, sep = ",")
colnames(lc_tl_z_mtx) <- rownames(lc_tl_z_mtx)
#cb_value <- min(-min(lc_tl_z_mtx), max(lc_tl_z_mtx))
png(paste0(out_path2, "lineage_zScore_heatmap.png"), width = 1100, height = 1100)
pheatmap(lc_tl_z_mtx, cluster_rows = T, cluster_cols = T, display_numbers = T, breaks = seq(-1,1,length.out = 100))
dev.off()
lc_tl_zscore_long <- export_for_cytoscape(lc_tl_z_mtx)
write.table(lc_tl_zscore_long, paste0(out_path2, "lc_tl_zScore_matrix_cytoscape.csv"), sep = ",", quote = F, row.names = F)


## QUANTIFICATION OF CLONAL DISTR ##
clone_id_pred_df <- data.frame(
   "cellID" = ferret_seurat$cellid,
   "cloneID" = ferret_seurat$cloneID,
   "celltype" = ferret_seurat$predicted.id
)
clone_id_pred_df <- clone_id_pred_df[!is.na(clone_id_pred_df$cloneID), ]


ct_per_clone <- sapply(unique(clone_id_pred_df$cloneID), function(cid) {clone_id_pred_df$celltype[clone_id_pred_df$cloneID == cid]})
num_cells_per_clone <- sapply(ct_per_clone, length)

table(num_cells_per_clone)
## ok
sum(table(num_cells_per_clone))
##
# Out of 9267 cells with recovered clone-IDs, 6349 cells were alone in their respective clone, i.e. the other cell of the clone were not recovered..
# 9267 - 6349 = 2918 cells in multi-cell clones
# 1038 multi-cell clones (1038/9267 = 0.1120104)
# 2918/1038 = 2.811175 cells per clone on average
# z-scores indicate the number of shared clones in respect to a random distribution;
# so pos z-score between RGCs and ENs indicate more-than random (i.e. biological) significance
##

## number of clones per trajectory:
traj1_cell_types <- c("(23) RGC S/G2/M","(11) RGC S/G2/M","(8) RGC G1","(18) RGC-IPC",
                      "(10) IPC G1","(2) IPC-Newborn EN","(19) IPC S/G2/M-Newborn EN","(3) IPC-Newborn EN","(0) IPC-Newborn EN",
                      "(7) General EN","(1) Upper-layer EN")
traj2_cell_types <- c("(22) RGC S/G2/M","(9) RGC S/G2/M","(17) RGC S/G2/M","(16) RGC G1","(4) RGC G1","(8) RGC G1","(18) RGC-IPC",
                      "(10) IPC G1","(2) IPC-Newborn EN","(19) IPC S/G2/M-Newborn EN","(3) IPC-Newborn EN","(0) IPC-Newborn EN",
                      "(7) General EN","(1) Upper-layer EN")
traj3_cell_types <- c("(22) RGC S/G2/M","(9) RGC S/G2/M","(14) IPC S/G2/M","(13) IPC S/G2/M","(24) IPC S/G2/M",
                      "(10) IPC G1","(2) IPC-Newborn EN","(19) IPC S/G2/M-Newborn EN","(3) IPC-Newborn EN","(0) IPC-Newborn EN",
                      "(7) General EN","(1) Upper-layer EN")

## only check multi-cell clones:
ct_per_clone_mcc <- ct_per_clone[num_cells_per_clone > 1]

traj1_count_ex <- 0; traj2_count_ex <- 0; traj3_count_ex <- 0
traj1_count_in <- 0; traj2_count_in <- 0; traj3_count_in <- 0
for(ct_vec in ct_per_clone_mcc) {
   in_traj_1 <- sum(ct_vec %in% traj1_cell_types)
   in_traj_2 <- sum(ct_vec %in% traj2_cell_types)
   in_traj_3 <- sum(ct_vec %in% traj3_cell_types)
   
   ## t1
   if(in_traj_1 == length(ct_vec)) {traj1_count_ex <- traj1_count_ex + 1}
   if(in_traj_1 > 0) {traj1_count_in <- traj1_count_in + 1}
   ## t2
   if(in_traj_2 == length(ct_vec)) {traj2_count_ex <- traj2_count_ex + 1}
   if(in_traj_2 > 0) {traj2_count_in <- traj2_count_in + 1}
   ## t1
   if(in_traj_3 == length(ct_vec)) {traj3_count_ex <- traj3_count_ex + 1}
   if(in_traj_3 > 0) {traj3_count_in <- traj3_count_in + 1}
}

traj1_count_ex; traj2_count_ex; traj3_count_ex
# [1] 479
# [1] 760
# [1] 511
traj1_count_in; traj2_count_in; traj3_count_in
# [1] 975
# [1] 1016
# [1] 963


## SAME FOR ONLY MITOTIC ##
traj1_cell_types_M <- c("(23) RGC S/G2/M","(11) RGC S/G2/M","(8) RGC G1","(18) RGC-IPC",
                        "(10) IPC G1","(19) IPC S/G2/M-Newborn EN")
traj2_cell_types_M <- c("(22) RGC S/G2/M","(9) RGC S/G2/M","(17) RGC S/G2/M","(16) RGC G1","(4) RGC G1","(8) RGC G1","(18) RGC-IPC",
                        "(10) IPC G1","(19) IPC S/G2/M-Newborn EN")
traj3_cell_types_M <- c("(22) RGC S/G2/M","(9) RGC S/G2/M","(14) IPC S/G2/M","(13) IPC S/G2/M","(24) IPC S/G2/M",
                        "(10) IPC G1","(19) IPC S/G2/M-Newborn EN")

traj1_count_ex_M <- 0; traj2_count_ex_M <- 0; traj3_count_ex_M <- 0
for(ct_vec in ct_per_clone_mcc) {
   in_traj_1 <- sum(ct_vec %in% traj1_cell_types_M)
   in_traj_2 <- sum(ct_vec %in% traj2_cell_types_M)
   in_traj_3 <- sum(ct_vec %in% traj3_cell_types_M)
   
   ## t1
   if(in_traj_1 == length(ct_vec)) {traj1_count_ex_M <- traj1_count_ex_M + 1}
   ## t2
   if(in_traj_2 == length(ct_vec)) {traj2_count_ex_M <- traj2_count_ex_M + 1}
   ## t1
   if(in_traj_3 == length(ct_vec)) {traj3_count_ex_M <- traj3_count_ex_M + 1}
}

traj1_count_ex_M; traj2_count_ex_M; traj3_count_ex_M
# [1] 44
# [1] 99
# [1] 28

cells_in_traj1_M <- sum(ferret_seurat$predicted.id %in% traj1_cell_types_M)
cells_in_traj2_M <- sum(ferret_seurat$predicted.id %in% traj2_cell_types_M)
cells_in_traj3_M <- sum(ferret_seurat$predicted.id %in% traj3_cell_types_M)

traj1_count_ex_M/cells_in_traj1_M * 100; traj2_count_ex_M/cells_in_traj2_M * 100; traj3_count_ex_M/cells_in_traj3_M * 100


##################### CLONE SIZE #############################
ipc_cell_types <- c("(2) IPC-Newborn EN","(13) IPC S/G2/M","(19) IPC S/G2/M-Newborn EN","(0) IPC-Newborn EN","(3) IPC-Newborn EN",
                    "(14) IPC S/G2/M","(18) RGC-IPC","(10) IPC G1","(24) IPC S/G2/M")
en_cell_types <- c("(7) General EN","(1) Upper-layer EN")
rgc_M_cell_types <- c("(23) RGC S/G2/M","(11) RGC S/G2/M","(22) RGC S/G2/M","(9) RGC S/G2/M","(17) RGC S/G2/M")
clone_in_ipc <- sapply(ct_per_clone_mcc, function(ct_vec) {
   if(sum(ct_vec %in% c(en_cell_types, rgc_M_cell_types)) == length(ct_vec)) {
      if(sum(ct_vec %in% en_cell_types) > 0 & sum(ct_vec %in% rgc_M_cell_types) > 0) {"direct"}
      else {"not_direct"}
   } else {"not_direct"}
})
sum(clone_in_ipc == "direct")
# [1] 113
sum(clone_in_ipc != "direct")
# [1] 925

clone_in_rgcA <- sapply(ct_per_clone_mcc, function(ct_vec) {
   if(sum(ct_vec %in% c(en_cell_types, c("(23) RGC S/G2/M","(11) RGC S/G2/M"))) == length(ct_vec)) {
      if(sum(ct_vec %in% en_cell_types) > 0 & sum(ct_vec %in% c("(23) RGC S/G2/M","(11) RGC S/G2/M")) > 0) {"direct"}
      else {"not_direct"}
   } else {"not_direct"}
})
sum(clone_in_rgcA == "direct")
# 28
## clones containing RGC_A:
sum(sapply(ct_per_clone_mcc, function(ct_vec) {sum(ct_vec %in% c("(23) RGC S/G2/M","(11) RGC S/G2/M")) > 1}))
# 38
## number of clones containing any RGC: 650

clone_in_rgcB <- sapply(ct_per_clone_mcc, function(ct_vec) {
   if(sum(ct_vec %in% c(en_cell_types, c("(22) RGC S/G2/M","(9) RGC S/G2/M"))) == length(ct_vec)) {
      if(sum(ct_vec %in% en_cell_types) > 0 & sum(ct_vec %in% c("(22) RGC S/G2/M","(9) RGC S/G2/M")) > 0) {"direct"}
      else {"not_direct"}
   } else {"not_direct"}
})
sum(clone_in_rgcB == "direct")
## 57


## plot cluster size:
num_cells_per_mcc <- sapply(ct_per_clone_mcc, length)
clone_size_plot_df <- data.frame(
   "number_cells_per_clone" = num_cells_per_mcc,
   "mode" = clone_in_ipc
)

library(ggsignif)
ggplot(clone_size_plot_df, aes(mode, number_cells_per_clone)) +
   geom_violin() +
   theme_bw() +
   geom_signif(comparisons = list(c("direct","not_direct")), map_signif_level = TRUE, test = "wilcox.test") +
   ylab("number of cells per clone")
ggsave("/datastore_share/Users/neuhaus/ferret_proj/results/pdf_plots/clone_size_direct_violin.pdf", width = 8, height = 6)

ggplot(clone_size_plot_df, aes(mode, log2(number_cells_per_clone))) +
   geom_violin() +
   theme_bw() +
   geom_signif(comparisons = list(c("direct","not_direct")), map_signif_level = TRUE, test = "wilcox.test") +
   ylab("log2(number of cells per clone)")
ggsave("/datastore_share/Users/neuhaus/ferret_proj/results/pdf_plots/clone_size_log2_direct_violin.pdf", width = 8, height = 6)

mean(clone_size_plot_df$number_cells_per_clone[clone_size_plot_df$mode == "direct"])
# [1] 2.389381
mean(clone_size_plot_df$number_cells_per_clone[clone_size_plot_df$mode == "not_direct"])
# [1] 2.862703

wilcox.test(number_cells_per_clone ~ mode, data = clone_size_plot_df)

