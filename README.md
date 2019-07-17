# SCC
SCC (scRNA-seq complementation). SCC can recover the gene expression of scRNA-seq data. There are two functions in SCC.

##SCC_filtration(data,celltype)
The data is the input matrix of scRNA-seq. The rownames represent the gene names and column represent cell names. The celltype represent the type of cells. This function can filtrate the outliers of scRNA-seq.

##SCC(data)
The only input data is the input matrix. The rownames represent the gene names and column represent cell names. SCC can recover the gene expression of scRNA-seq by a mixture model.
