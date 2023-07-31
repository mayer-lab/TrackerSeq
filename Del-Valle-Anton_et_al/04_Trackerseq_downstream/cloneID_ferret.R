## infer clone IDs from sparse matrix ##

library(readr)
#library(philentropy)
#library(dendextend)
library(data.table)
library(dplyr)
library(tibble)
library(igraph)
library(ggplot2)
library(ggrepel)

## function to create network from sparse matrix ##
build_graph_from_sparse_mtx <- function(mtx) {
  ## init isolated graph:
  g <- graph.empty(n = nrow(mtx), directed = F)
  V(g)$cellID <- rownames(mtx)
  ## save edges in matrix:
  ncells_per_bc <- apply(mtx, 2, sum)
  num_edges <- sapply(ncells_per_bc, function(n) {((n-1)*n)/2})
  edge_df <- matrix(nrow = sum(num_edges), ncol = 2); edge_attr_vec <- rep(NA, sum(num_edges))
  ## go through all lineage barcodes (columns) and save edges if cells are connected:
  i <- 1
  for(lbc_idx in colnames(mtx)) {
    lbc_idx <- as.character(lbc_idx)
    num_cells <- sum(mtx[,lbc_idx])
    if(num_cells == 0) {
      next
    } else if(num_cells == 1) {
      next
    } else if(num_cells == 2) {
      ## edge between two cells:
      g_idx <- which(mtx[, lbc_idx] == 1)
      edge_df[i,] <- c(g_idx[1], g_idx[2])
      edge_attr_vec[i] <- lbc_idx
      i <- i+1
    } else {
      ## num_cells >=3:
      g_idx <- which(mtx[, lbc_idx] == 1)
      for(j in 1:(length(g_idx) - 1)) {
        for(k in (j+1):length(g_idx)) {
          edge_df[i,] <- c(g_idx[j], g_idx[k])
          edge_attr_vec[i] <- lbc_idx
          i <- i+1
        }
      }
    }
  }
  edge_attr_vec <- edge_attr_vec[!is.na(edge_attr_vec)]
  edge_df <- edge_df[1:length(edge_attr_vec), ]
  
  edge_vec <- c(t(edge_df))
  
  ## update graph:
  g <- add_edges(g,
                 edges = edge_vec,
                 lbc_idx = edge_attr_vec
  )
  return(g)
}

## define output path:
out_path1 <- "/data/mayerlab/neuhaus/ferret_proj/final_results/00"
out_path2 <- "/datastore_share/Users/neuhaus/ferret_proj/final_results/00"

## read and prepare sparse matrix ##
sparse_matrix1 <- readr::read_csv("/data/mayerlab/neuhaus/ferret_proj/danger_zone/corrected_output_nread10_numi9/ED221124_sparse_matrix_UMI3_leven_read10_umi9.csv")
dim(sparse_matrix1)
sparse_matrix1$...1 <- NULL

## check integrity of matrix:
sum(is.na(sparse_matrix1$cellbc))
sum(is.null(sparse_matrix1$cellbc))

## not all cell-ids are unique (as we split them across 4 lanes)
length(unique(sparse_matrix1$cellbc))

## get rid of non-unique cell-bcs:
cell_bc_table1 <- table(sparse_matrix1$cellbc)
unique_cell_bcs <- names(cell_bc_table1)[as.numeric(cell_bc_table1) == 1]
unique_cell_bc_idx <- sparse_matrix1$cellbc %in% unique_cell_bcs
sparse_matrix1 <- sparse_matrix1[unique_cell_bc_idx, ]
length(unique(sparse_matrix1$cellbc)); dim(sparse_matrix1)

sparse_matrix1 <- column_to_rownames(sparse_matrix1, var = "cellbc")

## delete lineage barcodes that only existed in non-unique cells:
ncells_per_bc1 <- apply(sparse_matrix1, 2, sum)
length(ncells_per_bc1)
sparse_matrix1 <- sparse_matrix1[, ncells_per_bc1 > 0]

## look at sparse matrix:
sparse_matrix1[1:5,1:40]


## check the length of lineage barcodes:
barcode_list1 <- read.table("/data/mayerlab/neuhaus/ferret_proj/danger_zone/corrected_output_nread10_numi9/ED221124_barcode_list_leven_read10_umi9.txt")
head(barcode_list1)
barcode_vec1 <- barcode_list1$V1
## also subset for unique cell bcs:
barcode_vec1 <- barcode_vec1[unique_cell_bc_idx]
barcode_split1 <- lapply(barcode_vec1, function(el) {strsplit(el, ",")[[1]]})
head(barcode_split1)

barcode_single_vec1 <- unlist(barcode_split1)
## all lineage barcodes have length 37 (as they are supposed to):
table(sapply(barcode_single_vec1, nchar))


## check distributions of bc per cell:
ncells_per_bc <- apply(sparse_matrix1, 2, sum)
table(ncells_per_bc)
ncells_per_barcode_df <- data.frame(
  "cells_per_bc" = names(table(ncells_per_bc)),
  "count" = as.numeric(table(ncells_per_bc))
)
ggplot(ncells_per_barcode_df, aes(cells_per_bc, count, label = count)) +
  geom_point() +
  geom_text_repel() +
  theme_bw() +
  xlab("cells per lineage barcode")
ggsave(paste0(out_path2, "/ncells_per_lbc_scatter.png"), width = 10, height = 8)
ggsave(paste0(out_path1, "/ncells_per_lbc_scatter.png"), width = 10, height = 8)

## number of lineage barcodes per cell:
nlbc_per_cell <- apply(sparse_matrix1, 1, sum)
table(nlbc_per_cell)
nlbc_per_cell_df <- data.frame(
  "lbc_per_cell" = as.numeric(names(table(nlbc_per_cell))),
  "count" = as.numeric(table(nlbc_per_cell))
)
nlbc_per_cell_df$lbc_per_cell <- factor(nlbc_per_cell_df$lbc_per_cell, levels = as.numeric(nlbc_per_cell_df$lbc_per_cell))
ggplot(nlbc_per_cell_df, aes(lbc_per_cell, count, label = count)) +
  geom_point() +
  geom_text_repel() +
  theme_bw() +
  xlab("lineage barcodes per cell")
ggsave(paste0(out_path2, "/nlbc_per_cell_scatter.png"), width = 10, height = 8)
ggsave(paste0(out_path1, "/nlbc_per_cell_scatter.png"), width = 10, height = 8)



## build graph from sparse matrix ##
g <- build_graph_from_sparse_mtx(mtx = sparse_matrix1)
saveRDS(g, "/data/mayerlab/neuhaus/ferret_proj/final_results/00/cloneID_igraph_nread10_numi9.rds")

## calculate connected components ##
g_comp <- components(g)

cells_per_clone <- table(g_comp$membership)
table(cells_per_clone)

ncells_per_clone_df <- data.frame(
  "cells_per_clone" = as.numeric(names(table(cells_per_clone))),
  "count" = as.numeric(table(cells_per_clone))
)
ggplot(ncells_per_clone_df, aes(cells_per_clone, count, label = count)) +
  geom_point() +
  geom_text_repel() +
  theme_bw() +
  xlab("cells per clone")
ggsave(paste0(out_path2, "/ncells_per_clone_scatter.png"), width = 10, height = 8)
ggsave(paste0(out_path1, "/ncells_per_clone_scatter.png"), width = 10, height = 8)

#table(degree(g))

clones <- data.frame(
  "cloneID" = g_comp$membership,
  "cellbc" = V(g)$cellID
)

## remove lineage information of large clone
huge_clone_id <- which(g_comp$csize == 4175)
clones$cloneID_filtered <- clones$cloneID
clones$cloneID_filtered[clones$cloneID_filtered == huge_clone_id] <- NA

## write clonal information:
write.csv(clones, file = paste0(out_path1, "/ED221124_cloneIDs_nreads10_numi9.csv"))

