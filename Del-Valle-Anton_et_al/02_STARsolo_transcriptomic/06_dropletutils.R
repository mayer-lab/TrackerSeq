library(tidyverse)
library(DropletUtils)
library(Seurat)

dir.name <- list("/usr/users/mmnb1ngs/collab_work/workspace/ykbabal/Novaseq_Ferret/res_STAR/T1/Solo.out/Gene/raw",
                 "/usr/users/mmnb1ngs/collab_work/workspace/ykbabal/Novaseq_Ferret/res_STAR/T2/Solo.out/Gene/raw",
                 "/usr/users/mmnb1ngs/collab_work/workspace/ykbabal/Novaseq_Ferret/res_STAR/T3/Solo.out/Gene/raw",
                 "/usr/users/mmnb1ngs/collab_work/workspace/ykbabal/Novaseq_Ferret/res_STAR/T4/Solo.out/Gene/raw")

list.files(dir.name[[1]])

## Sample 1
sce1 <- read10xCounts(dir.name[[1]])
sce1

class(counts(sce1))

my.counts <- counts(sce1)


e.out <- emptyDrops(my.counts, niters=100000, lower = 500)
is.cell <- e.out$FDR <= 0.001
sum(is.cell, na.rm=TRUE)


sce1.filtered <- sce1[,which(is.cell)]

class(sce1.filtered)
colnames(sce1.filtered) <- sce1.filtered$Barcode

filtered.counts <- counts(sce1.filtered)
dim(filtered.counts)

write10xCounts(path = "/usr/users/mmnb1ngs/collab_work/workspace/ykbabal/Novaseq_Ferret/dropletutils_output/T1", 
               x = filtered.counts, type = "sparse")

suerat <- Read10X(data.dir = "/usr/users/mmnb1ngs/collab_work/workspace/ykbabal/Novaseq_Ferret/dropletutils_output/T1")
s1 <- CreateSeuratObject(suerat)
s1

saveRDS(s1, file = "/usr/users/mmnb1ngs/collab_work/workspace/ykbabal/Novaseq_Ferret/seurat_objects/T1_seurat.RDS")


## Sample 2
sce1 <- read10xCounts(dir.name[[2]])
sce1

class(counts(sce1))

my.counts <- counts(sce1)


e.out <- emptyDrops(my.counts, niters=100000, lower = 500)
is.cell <- e.out$FDR <= 0.001
sum(is.cell, na.rm=TRUE)


sce1.filtered <- sce1[,which(is.cell)]

class(sce1.filtered)
colnames(sce1.filtered) <- sce1.filtered$Barcode

filtered.counts <- counts(sce1.filtered)
dim(filtered.counts)

write10xCounts(path = "/usr/users/mmnb1ngs/collab_work/workspace/ykbabal/Novaseq_Ferret/dropletutils_output/T2", 
               x = filtered.counts, type = "sparse")

suerat <- Read10X(data.dir = "/usr/users/mmnb1ngs/collab_work/workspace/ykbabal/Novaseq_Ferret/dropletutils_output/T2")
s1 <- CreateSeuratObject(suerat)
s1

saveRDS(s1, file = "/usr/users/mmnb1ngs/collab_work/workspace/ykbabal/Novaseq_Ferret/seurat_objects/T2_seurat.RDS")


## Sample 3
sce1 <- read10xCounts(dir.name[[3]])
sce1

class(counts(sce1))

my.counts <- counts(sce1)


e.out <- emptyDrops(my.counts, niters=100000, lower = 500)
is.cell <- e.out$FDR <= 0.001
sum(is.cell, na.rm=TRUE)


sce1.filtered <- sce1[,which(is.cell)]

class(sce1.filtered)
colnames(sce1.filtered) <- sce1.filtered$Barcode

filtered.counts <- counts(sce1.filtered)
dim(filtered.counts)

write10xCounts(path = "/usr/users/mmnb1ngs/collab_work/workspace/ykbabal/Novaseq_Ferret/dropletutils_output/T3", 
               x = filtered.counts, type = "sparse")

suerat <- Read10X(data.dir = "/usr/users/mmnb1ngs/collab_work/workspace/ykbabal/Novaseq_Ferret/dropletutils_output/T3")
s1 <- CreateSeuratObject(suerat)
s1

saveRDS(s1, file = "/usr/users/mmnb1ngs/collab_work/workspace/ykbabal/Novaseq_Ferret/seurat_objects/T3_seurat.RDS")

## Sample 4
sce1 <- read10xCounts(dir.name[[4]])
sce1

class(counts(sce1))

my.counts <- counts(sce1)


e.out <- emptyDrops(my.counts, niters=100000, lower = 500)
is.cell <- e.out$FDR <= 0.001
sum(is.cell, na.rm=TRUE)


sce1.filtered <- sce1[,which(is.cell)]

class(sce1.filtered)
colnames(sce1.filtered) <- sce1.filtered$Barcode

filtered.counts <- counts(sce1.filtered)
dim(filtered.counts)

write10xCounts(path = "/usr/users/mmnb1ngs/collab_work/workspace/ykbabal/Novaseq_Ferret/dropletutils_output/T4", 
               x = filtered.counts, type = "sparse")

suerat <- Read10X(data.dir = "/usr/users/mmnb1ngs/collab_work/workspace/ykbabal/Novaseq_Ferret/dropletutils_output/T4")
s1 <- CreateSeuratObject(suerat)
s1

saveRDS(s1, file = "/usr/users/mmnb1ngs/collab_work/workspace/ykbabal/Novaseq_Ferret/seurat_objects/T4_seurat.RDS")



