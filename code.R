library(SingleCellExperiment)
library(scmap)
library(mclust)
##读入GSM2230757的数据
data<-read.csv("C:/Users/zy/Desktop/scRNAdata/GSM2230757/GSM2230757_human1_umifm_counts.csv")
rownames(data) <- data[,1]
data<-data[,-1]
data<-data[,-1]
celltype<-data[,1]
data<-data[,-1]
data<-t(data)
#scmapcell
sce <- SingleCellExperiment(assays = list(normcounts = as.matrix(data)), colData = celltype)
logcounts(sce) <- log2(normcounts(sce) + 1)

rowData(sce)$feature_symbol <- rownames(sce)
isSpike(sce, "ERCC") <- grepl("^ERCC-", rownames(sce))

sce <- sce[!duplicated(rownames(sce)), ]
sce <- selectFeatures(sce, suppress_plot = FALSE)
set.seed(1)
sce <- indexCell(sce)
scmapCell_results <- scmapCell(
  sce, 
  list(
    data = metadata(sce)$scmap_cell_index
  )
)
#合并数据
cell<-scmapCell_results$data$cells
sim<-scmapCell_results$data$similarities
data1<-data
celllen<-rep(0,ncol(data1))
for (i in 1:ncol(data1))
{
  print(i)
  for (j in 1:nrow(data1))
  {
    celllen[i]<-celllen[i]+data1[j,i]*data1[j,i]
  }
  celllen[i]<-sqrt(celllen[i])
  data1[,i]<-data1[,i]*100000/celllen[i]
}
ad<-matrix(nrow=nrow(data1),ncol=10)
for (i in 1:ncol(data1))
{
  nearpoint <- rep(NA,10)
  for (j in 1:10)
  {
    nearpoint[j]<-cell[j,i]
    ad[,j]<-data1[,nearpoint[j]]
  }
}