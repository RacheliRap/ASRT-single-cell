############################################################################################
# sci-RNA-seq analysis.
# Based on https://www.singlecellcourse.org/single-cell-rna-seq-analysis-using-seurat.html
# The data is from https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE184602
# The data was filtetred before published, therefore we will skip the filtering part
############################################################################################


library(Seurat)
library(ggplot2)
library(SingleR)
library(dplyr)
library(celldex)
library(RColorBrewer)
library(SingleCellExperiment)

setwd("C:/Users/rachelRa/")

file_names <- dir("Documents/sci-RNA-seq/allelic_sci-RNA-seq_catenated_count_matrices/", full.names = T) #where you have your files

adj.matrix <- do.call(rbind,lapply(file_names,read.delim))

#adj.matrix <- Read10X("Documents/sci-RNA-seq/allelic_sci-RNA-seq_catenated_count_matrices/")
adj.matrix <- data.frame(apply(adj.matrix[,-1], 2, as.numeric), row.names=adj.matrix[,1])

annotation = t(rbind.data.frame(strsplit(row.names(adj.matrix), "_")))
annotation.allele1 = annotation[grepl("ref", annotation[,1]),]
annotation.allele2 = annotation[grepl("alt", annotation[,1]),]

adj.matrix.allelle1 <- adj.matrix[grepl("ref", row.names(adj.matrix)), ]
adj.matrix.allelle2 <- adj.matrix[grepl("alt", row.names(adj.matrix)), ]

adj.matrix.allelle1.NPC = adj.matrix.allelle1[,grep("NPC", names(adj.matrix.allelle1))]
adj.matrix.allelle2.NPC = adj.matrix.allelle2[,grep("NPC", names(adj.matrix.allelle2))]


row.names(adj.matrix.allelle1.NPC) <- annotation.allele1[,2]
row.names(adj.matrix.allelle2.NPC) <- annotation.allele2[,2]

is.mono = NULL
propA = NULL
is.not.bi = NULL
THRS = 0
for (i in 1:nrow(adj.matrix.allelle1.NPC)) {
  geneA = adj.matrix.allelle1[i,]
  geneB = adj.matrix.allelle2[i,]
  is.mono[i] = sum(geneA>THRS & geneB>THRS) == 0 & sum(xor(geneA>THRS, geneB>THRS)) >= 20
  is.not.bi[i] = sum(geneA>THRS & geneB>THRS) == 0
  #propA = rbind.data.frame(propA, as.numeric(geneA>THRS) - as.numeric(geneB>THRS))
}

adj.matrix.allelle1.mono = adj.matrix.allelle1.NPC[is.not.bi,]



srat <- CreateSeuratObject(adj.matrix.allelle1.mono,project = "pbmc10k") 

adj.matrix <- NULL
str(srat)

meta <- srat@meta.data

summary(meta$nCount_RNA)
summary(meta$nFeature_RNA)

srat <- NormalizeData(srat)

srat <- FindVariableFeatures(srat, selection.method = "vst", nfeatures = 100)

top10 <- head(VariableFeatures(srat), 10)
top10

plot1 <- VariableFeaturePlot(srat)
LabelPoints(plot = plot1, points = top10, repel = TRUE, xnudge = 0, ynudge = 0)

all.genes <- rownames(srat)
srat <- ScaleData(srat, features = all.genes)

srat <- RunPCA(srat, features = VariableFeatures(object = srat))
ElbowPlot(srat)


srat <- FindNeighbors(srat, dims = 1:10)
srat <- FindClusters(srat, resolution = 0.5)

srat <- RunUMAP(srat, dims = 1:10, verbose = F)


DimPlot(srat,label.size = 4,repel = T,label = T)
