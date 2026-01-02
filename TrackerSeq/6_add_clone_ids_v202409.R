library(data.table)
library(dplyr)
library(tibble)
library(igraph)
library(tidyverse)
library(philentropy)
library(dendextend)
library(parallel)
library(Seurat)
library(ggplot2)
library(reshape2)
library(viridis)
library(rgl)
library(data.table)
library(waldo)
library(scCustomize)
library(prismatic)
library(RColorBrewer)
library(ggrepel)
library(MASS)
library(svglite)


# input cellbc-lbc matrix
matrix_csv<-'<input-path>'

########################################################################################
# Function for matching clone list cellbcs with seurat cellbcs and adding metadata column "cloneID"
merge.transcr.ling.may.pB2 <- function(seuratobject, pool, pool.col = "cloneID", seurat.new.col = "cloneID", ... ){
  cells <- intersect(pool$cellbc,colnames(seuratobject))
  rownames(pool) <- pool$cellbc
  seuratobject@meta.data[, seurat.new.col] <- "No LBC"
  seuratobject@meta.data[cells, seurat.new.col] <- pool[cells, pool.col]
  Idents(object = seuratobject) <- seurat.new.col
  seuratobject@meta.data[WhichCells(seuratobject, idents = "No LBC"), seurat.new.col] <- NA
  return(seuratobject)
}
########################################################################################
########################################################################################
## function to create network with pattern weighting information ##
## Note that this script is more computationally intensive, taking ~1-24hrs @ single core #
build_graph_from_sparse_mtx_w_pattern_weights <- function(mtx, n_cores) {
  ## init isolated graph:
  g <- graph.empty(n = nrow(mtx), directed = F)
  V(g)$cellID <- rownames(mtx)
  
  # iterate through non-zero barcodes
  clean_mtx<-mtx[, which(colSums(mtx)>1)]
  lbc_vec<-colnames(clean_mtx)
  
  # parrallelize the edge finding and weight calculation process
  edge_output <- mcmapply(lbc_vec, mc.cores = n_cores, FUN = function(lbc_idx){
    lbc_idx <- as.character(lbc_idx)
    g_idx <- which(mtx[, lbc_idx] == 1)
    g_idx_mtx<-mtx[g_idx,]
    ncells_per_bc <- nrow(g_idx_mtx)
    num_edges <- ((ncells_per_bc-1)*ncells_per_bc)/2
    edge_df <- matrix(nrow = sum(num_edges), ncol = 2); edge_attr_vec <- rep(NA, sum(num_edges)); weight_df <- rep(NA, sum(num_edges))
    #which(colnames(mtx) %in% mtx[g_idx,]),]
    i <- 1
    for(j in 1:(length(g_idx) - 1)) {
      # for each pairwise comparison of lbcs, if they match create an edge, if not move to next pairwise set
      # weight the edge by the % of barcodes contained within each cell
      for(k in (j+1):length(g_idx)) {
        bc_frac<-2*length(intersect(which(g_idx_mtx[j,] == 1), which(g_idx_mtx[k,] == 1)))/(as.numeric(length(unique(which(g_idx_mtx[j,] == 1)))) + as.numeric(length(unique(which(g_idx_mtx[k,] == 1)))))
        edge_df[i,] <- c(g_idx[j],  g_idx[k])
        edge_attr_vec[i] <- lbc_idx
        weight_df[i] <- bc_frac
        i <- i + 1
      }
    }
    outputter<-list(edge_df, edge_attr_vec, weight_df)
    return(outputter)
  }
  )
  edge_output <- t(as.data.frame(edge_output))
  edge_df <- do.call(rbind, edge_output[,1])
  edge_attr_vec <- unlist(edge_output[,2])
  weight_df <- unlist(edge_output[,3])
  edge_vec <- c(t(edge_df))
  weight_vec <-c(t(weight_df))
  
  g <- add_edges(g,
                 edges = edge_vec,
                 lbc_idx = edge_attr_vec,
                 weight = weight_vec
  )
  return(g)
}

########################################################################################

########################################################################################
# Apply threshold to obtain filtered network with patterns above threshold# 
# Note: threshold is inclusive as implemented
network_patterning <- function(gr, threshold){
  if (threshold == 1){
    exclusion_edges <- E(gr)[weight < threshold]
  }
  else if (threshold < 1){
    exclusion_edges <- E(gr)[weight <= threshold]
    
  }
  filtered_graph <- delete_edges(gr, exclusion_edges)
  filtered_components <- igraph::components(filtered_graph)
  clones <- na.omit(data.frame(
    "cloneID" = filtered_components$membership,
    "cellbc" = V(filtered_graph)$cellID
  ))
  return(clones)
}
########################################################################################


########################################################################################
# Iterate over many thresholds (set by test_thresholds vector) to determine the shoulder in the clone ID based on the weight threshold!? 
iterate_clone_pattern <- function(gr, test_thresholds){
  # initialize with perfect pattern matches only
  h <- 2
  all_clones <- network_patterning(gr, threshold = 1)
  clone_names <- c("cellbc", "cloneID_1")
  # iterate through the test_thresholds list
  for (i in test_thresholds){
    h <- h + 1
    clones <- network_patterning(gr, i)
    clone_names[h] <- paste("cloneID", i, sep = "_")
    all_clones <- merge(all_clones, clones, by = "cellbc", all = T, 
                        suffixes = c(clone_names[h-1], clone_names[h]))
  }
  colnames(all_clones) <- clone_names
  return(all_clones)
}
########################################################################################

# _______MAIN WORKING CODE SECTION _____ # 
## Import cellbc-lbc matrix and lbc barcodes files, compute distribution statistics for filtering out high-lbc cells, 
## which may create bundling errors

# Import the sparse clone matrix from python pipeline
sparse_matrix <- readr::read_csv(file = matrix_csv, col_names = T)
sparse_matrix <- as.data.frame(sparse_matrix)
rownames(sparse_matrix)<-sparse_matrix$cellbc
sparse_matrix$cellbc <- NULL

# compute total lbcs per cell and cells per lbc (col and row sums)
sparse_matrix_rowsums<-rowSums(sparse_matrix)
sparse_matrix_rowsums<-as.data.frame(sparse_matrix_rowsums)
sparse_matrix_colsums<-colSums(sparse_matrix)
sparse_matrix_colsums<-as.data.frame(sparse_matrix_colsums)

statz<-quantile(sparse_matrix_rowsums$sparse_matrix_rowsums, probs = c(0.85,0.9,0.95))

# Fit the distribution of lbcs per cell using exponential regression
mean(sparse_matrix_rowsums$sparse_matrix_rowsums)
model_dist<- MASS::fitdistr(sparse_matrix_rowsums$sparse_matrix_rowsums,"Exponential")

lambda <- 1/model_dist$estimate
inverse_lambda<- model_dist$estimate
# Set functions for dist model results
pois_funct<-function(x){length(sparse_matrix_rowsums$sparse_matrix_rowsums)*((lambda^x*(exp(-lambda)))/factorial(x))}
exp_funct<-function(x){length(sparse_matrix_rowsums$sparse_matrix_rowsums)*(inverse_lambda*(exp(-1*inverse_lambda*x)))}

# Plot hist of lbcs per cell and stats
ggplot(sparse_matrix_rowsums, aes(x=sparse_matrix_rowsums)) + 
  geom_histogram(binwidth=1, fill="black") +
  stat_function(fun = exp_funct, color = "magenta", alpha = 0.7) +
  xlab("# LBCs per cell") + 
  ylab("# cells containing") +
  ggtitle("Modeling LBCs per cell: exponential fit") #+


ggplot(sparse_matrix_colsums, aes(x=sparse_matrix_colsums)) + 
  geom_histogram(binwidth=1, fill="black") +
  xlab("# cells per LBC") + 
  ylab("# lbcs") +
  ggtitle("Modeling cells per LBC")


# Generate network (igraph object) from the sparse matrix with edges for each lbc shared bewteen two cells (vertices). 
# Each edge is weighted by a calculation of the proportion of shared barcodes bewteen the two cells
# Process runs multicore - can still take quite awhile especially for large cell-lbc matrices (e.g. >500x500)
g <- build_graph_from_sparse_mtx_w_pattern_weights(sparse_matrix, n_cores = 16)
# recommend saving global environment variables at this stage!

# Iteratively filter the graph using pattern match thresholds (e.g. 1000 cutoffs between 0 and 1, can be adjusted)
# Note the fractional calculations typically only have certain values which creates a stepping effect in resulting clone size lists / statistics...
test_thresholds <- seq(from = 0, to = 0.999 , by = 0.001)
#iterative_clone_df <- iterate_clone_pattern(g, test_thresholds)
iterative_clone_df <- iterate_clone_pattern(g, test_thresholds)
iterative_lists <- colnames(iterative_clone_df)
iterative_lists<-as.data.frame(iterative_lists[2:length(iterative_lists)])
colnames(iterative_lists)<-c("run")

# Generate descriptive data frame for all clone lists (which are otherwise not exported from above code)
# should convert into an apply function for efficiency... but it's not too bad atm
for (i in 1:length(iterative_lists$run)){
  #print(i)
  l <- length(unique(iterative_clone_df[, i+1]))
  iterative_lists$length[i]<-l
  data_name<-iterative_lists$run[i]
  x<- iterative_clone_df %>% dplyr::count(iterative_clone_df[i+1])
  x_size <- max(x$n)
  multicells<-filter(x, n>1)
  t5_threshold<-quantile(multicells$n, probs = c(0.95))
  t10_threshold<-quantile(multicells$n, probs = c(0.90))
  t25_threshold<-quantile(multicells$n, probs = c(0.75))
  top5_percentile<-filter(x, n>t5_threshold)
  top10_percentile<-filter(x, n>t10_threshold)
  top25_percentile<-filter(x, n>t25_threshold)
  singles <- filter(x, n==1)
  n_single <- length(singles$n)
  n_multi <-length(multicells$n)
  x_avg <- mean(multicells$n)
  top5_avg<-mean(top5_percentile$n)
  top10_avg<-mean(top10_percentile$n)
  top25_avg<-mean(top25_percentile$n)
  x_median <-median(multicells$n)
  x_variance <-var(multicells$n)
  top5_variance<-var(top5_percentile$n)
  top10_variance<-var(top10_percentile$n)
  top25_variance<-var(top25_percentile$n)
  iterative_lists$biggest_clone_size[i]<-x_size
  iterative_lists$number_single[i]<-n_single
  iterative_lists$number_multicells[i]<-n_multi
  iterative_lists$avg_clone_size[i]<-x_avg
  iterative_lists$top5_avg[i]<-top5_avg
  iterative_lists$top10_avg[i]<-top10_avg
  iterative_lists$top25_avg[i]<-top25_avg
  iterative_lists$median_clone_size[i]<-x_median
  iterative_lists$clone_size_variance[i]<-x_variance
  iterative_lists$top5_variance[i]<-top5_variance
  iterative_lists$top10_variance[i]<-top10_variance
  iterative_lists$top25_variance[i]<-top25_variance
}

# Data massage
iterative_lists$threshold<-parse_number(iterative_lists$run)
iterative_lists$biggest_clone_size<-as.numeric(iterative_lists$biggest_clone_size)
iterative_lists$avg_clone_size<-as.numeric(iterative_lists$avg_clone_size)
iterative_lists$top5_avg<-as.numeric(iterative_lists$top5_avg)
iterative_lists$top10_avg<-as.numeric(iterative_lists$top10_avg)
iterative_lists$top25_avg<-as.numeric(iterative_lists$top25_avg)
iterative_lists$median_clone_size<-as.numeric(iterative_lists$median_clone_size)
iterative_lists$clone_size_variance<-as.numeric(iterative_lists$clone_size_variance)
iterative_lists$top5_variance<-as.numeric(iterative_lists$top5_variance)
iterative_lists$top10_variance<-as.numeric(iterative_lists$top10_variance)
iterative_lists$top25_variance<-as.numeric(iterative_lists$top25_variance)

# Save the iterative clone df for future reference
write.csv(iterative_clone_df, '<filepath>')



# Pattern statistics visualization
weight_df <- E(g)$weight
weight_df <- as.data.frame(weight_df)
ggplot(weight_df, aes(x=weight_df)) + 
  geom_histogram(binwidth=0.005, fill = "black") +
  theme(plot.title = element_text(size=11, face = 'bold', family = "Arial")) +
  ggtitle("Histogram of Pattern scores") +
  xlab("Pattern match threshold") +
  ylab("Number of cell pairs")

ggsave(device = "png", path = "<directory>", dpi = 1000, filename = "<filename>.png")

# examine clones based on inclusion of all network edges, regardless of weight
unpatterned<-iterative_clone_df %>% dplyr::count(iterative_clone_df["cloneID_0"])
unpatterned<-filter(unpatterned, n>1)
unpatterned <- unpatterned[order(unpatterned$n, decreasing = T),]
rownames(unpatterned)<-1:nrow(unpatterned)
ggplot(unpatterned, aes(x=as.numeric(rownames(unpatterned)), y=n)) +
  geom_area(color = 1, lwd = 1) +
  ylab("Cells per clone") +
  xlab("Unique clone") +
  ggtitle("Clone sizes resulting from direct LBC match")



# Plot clone sizes as the pattern match threshold is adjusted from 0 to 1
ggplot(iterative_lists) + 
  #geom_point(aes(x = threshold, y = length), size = 0.1, color = "black") +
  #geom_text(x=0.75, y=2000, label="# unique clones") +
  geom_point(aes(x = threshold, y = number_multicells), size = 0.5, color = "steelblue4") +
  geom_text(x=0.75, y=60, label="# unique multicell clones", color = "steelblue4") +
  geom_point(aes(x = threshold, y = biggest_clone_size), size = 0.5, color = "red", fill = "red") +
  geom_text(x=0.75, y=45, label="size of largest clone", color = "red") +
  #geom_point(aes(x = threshold, y = number_single), size = 0.1, color = "chartreuse4") +
  #geom_text(x=0.75, y=1800, label="# unique singlecell clones", color = "chartreuse4") +
  theme(plot.title = element_text(size=11, face = 'bold', family = "Arial")) +
  xlab("Pattern match threshold") +
  ylab('Number of clones') +
  ggtitle('Clone dynamics across pattern assignment spectrum')


# Zoom in on largest clones
ggplot(iterative_lists) + 
  geom_point(aes(x = threshold, y = biggest_clone_size), size = 0.1, color = "red") +
  geom_text(x=0.75, y = 30, label="size of largest clone", color = "red") +
  xlab("Pattern match threshold") +
  ylab('Number of cells in clone') +
  ggtitle('Clone size dynamics across pattern assignment spectrum') +
  theme(plot.title = element_text(size=11, face = 'bold', family = "Arial")) +
  #geom_point(aes(x = threshold, y = biggest_clone_size))
  coord_cartesian(xlim=NULL, ylim=c(0, 50))

# Plot multicell clone stats
ggplot(iterative_lists) + 
  geom_point(aes(x = threshold, y = avg_clone_size), size = 0.5, color = "blue") +
  geom_text(x=0.75, y=8, label="mean multicell clone size", color = "blue") +
  #geom_point(aes(x = threshold, y = top5_avg), size = 0.1, color = "firebrick3") +
  #geom_text(x=0.75, y=25, label="95th percentile mean multicell clone size", color = "firebrick3") +
  geom_point(aes(x = threshold, y = top10_avg), size = 0.5, color = "purple") +
  geom_text(x=0.75, y=21, label="90th percentile mean multicell clone size", color = "purple") +
  #geom_point(aes(x = threshold, y = top25_avg), size = 0.1, color = "darkcyan") +
  #geom_text(x=0.75, y=21, label="75th percentile mean multicell clone size", color = "darkcyan") +
  geom_point(aes(x = threshold, y = median_clone_size), size = 0.5, color = "darkgreen") +
  geom_text(x=0.75, y=13, label="median multicell clone size", color = "darkgreen") +
  xlab("Pattern match threshold") +
  ylab('Value') +
  ggtitle('Mean & median statistics across pattern assignment spectrum') +
  theme(plot.title = element_text(size=11, face = 'bold', family = "Arial")) +
  coord_cartesian(xlim=NULL, ylim=c(0, 25))

# Plot standard deviations
ggplot(iterative_lists) + 
  geom_point(aes(x = threshold, y = clone_size_variance^.5), size = 0.5, color = "black") +
  geom_text(x=0.5, y=15, label="all multicell clones", color = "black") +
  #geom_point(aes(x = threshold, y = top25_variance^.5), size = 0.5, color = "blue") +
  #geom_text(x=0.5, y=17.5, label="75th percentile multicell clones", color = "blue") +
  #geom_point(aes(x = threshold, y = top10_variance^.5), size = 0.5, color = "darkgreen") +
  #geom_text(x=0.5, y=20, label="90th percentile multicell clones", color = "darkgreen") +
  #geom_point(aes(x = threshold, y = top5_variance^.5), size = 0.5, color = "darkorange3") +
  #geom_text(x=0.5, y=22.5, label="95th percentile multicell clones", color = "darkorange3") +
  theme(plot.title = element_text(size=11, face = 'bold', family = "Arial")) +
  xlab("Pattern match threshold") +
  ylab('Standard deviation') +
  ggtitle('Standard deviation of multicell clones across pattern assignment spectrum') #+
  #coord_cartesian(xlim=NULL, ylim=c(0, 60))

# Zoom in on st. dev statistic
ggplot(iterative_lists) + 
  geom_point(aes(x = threshold, y = clone_size_variance), size = 0.5, color = "black") +
  xlab("Pattern match threshold") +
  ylab('Variance') +
  ggtitle('Variance of multicell clones across pattern assignment spectrum') +
  theme(plot.title = element_text(size=11, face = 'bold', family = "Arial")) #+
  #coord_cartesian(xlim=NULL, ylim=c(0, 10))


# Based on plots generated in the above section, select a pattern threshold that optimizes for clone resolution and reduciton of noise
# After selection of the desired threshold for LBC pattern cutoff, generate a single clone list. 

# name the clondID column that corresponds to the desired cutoff level:
sparse_matrix_clones<-iterative_clone_df[, c('cellbc', 'cloneID_0.445')]

sparse_matrix_dt<-as.data.table(sparse_matrix_clones)
colnames(sparse_matrix_dt)<-c('cellbc', 'cloneID')

# compute the # of cells in each clone
sparse_matrix_clone_count<-sparse_matrix_dt[, .N, by = .(cloneID)]
# filter for multicells
sparse_matrix_clone_count<-filter(sparse_matrix_clone_count, cloneID>1)

# save the clone list
write.csv(sparse_matrix_clones, '<filepath>')

# prepare to match into Seurat object, which may involve string edits to the cellbcs
# copy clone list before making edits to cellbcs
sparse_matrix_clones_output<- sparse_matrix_clones
colnames(sparse_matrix_clones_output)<-c('cellbc', 'cloneID')
sparse_matrix_clones_output$sample<-"<samplename>"
# Example edits to cloneIDs (in the case of multiple TrackerSeq libraries per Seurat object, to enforce distinct clones)
sparse_matrix_clones_output$cloneID <- paste0(sparse_matrix_clones_output$cloneID, "<suffix>")
# Example edits to cellbcs to match those in Seurat (which usually have some suffix, potentially depending on sample in cases of samples with same cellbc from different GEM wells)
sparse_matrix_clones_output$cellbc <- paste("<prefix>", sparse_matrix_clones_output$cellbc, sep = "")

# If multiple cell - lbc lineage matrices, must merge clone lists prior to merging to merged seurat object
merged_pool<-do.call("rbind", list(sparse_matrix_clones_2, sparse_matrix_clones_3, sparse_matrix_clones_4, sparse_matrix_clones_5))
write.csv(merged_pool, '<filepath>')

# Read in Seurat
seurat_object <- readRDS('<seurat_filepath>')

# Match-in clone file to cloneID field in Seurat
seurat_object@meta.data$cloneID <- NA
seurat_object <- merge.transcr.ling.may.pB2(seurat_object, pool = merged_pool, pool.col = "cloneID", seurat.new.col = "cloneID")
seurat_object@meta.data$cloneID<-seurat_object@meta.data$cloneID

# Pull seurat metadata into data table for quick counting using data table
dt <- seurat_object@meta.data %>% as.data.table
# Clonelist matched to seurat, multicell clone list
seurat_clones <- na.omit(dt[, .N, by = .(cloneID, TSlib)])
multicell_clones <- na.omit(seurat_clones[seurat_clones$N>1])

seurat_clones_bycluster <- na.omit(dt[, .N, by = .(cloneID, TSlib, seurat_clusters)])
multicell_seurat_clones_bycluster <- na.omit(seurat_clones_bycluster[seurat_clones_bycluster$N>1])
cluster_multicell_counts<-na.omit(multicell_seurat_clones_bycluster[, .N, by = .(seurat_clusters)])

# Clone histogram
multicell_clones<- arrange(multicell_clones, desc(N))
ggplot(multicell_clones, aes(x=parse_number(rownames(multicell_clones)), y = N)) + 
  geom_col(just = 0, width = 1, fill = "black") + 
  xlim(0, 150) + 
  ylim(0, 20) + 
  xlab("Unique clones") + 
  ylab("Number of cells") +
  ggtitle("Multicell clones, matched into seurat")


# Add clone category field in seurat, plot multicell clones in UMAP
# can set the non-barcoded cells to NA if desired for plotting, but often useful to have a character value here for downstream subsetting
seurat_object@meta.data$multicell_clones<-"No LBC"
seurat_object@meta.data$multicell_clones[is.na(seurat_object@meta.data$cloneID)==FALSE]<-"Singlecell"
seurat_object@meta.data$multicell_clones[seurat_object@meta.data$cloneID%in%e246_multicell_clones$cloneID]<-"Multicell"
DimPlot(seurat_object, reduction = "umap", label = FALSE, repel = TRUE, group.by = "multicell_clones") + scale_color_manual(values=c("red", "grey50", "blue")) +ggtitle("") + labs(x=NULL, y=NULL)

# Save out the seurat file with it's new cloneIDs! 
saveRDS(seurat_object, seurat_RDS)






