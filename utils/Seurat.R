library(Seurat)
library(reticulate)

use_python("/home/bruce1996/anaconda3/bin/python")
pd <- import("pandas")
df <- pd$read_pickle("/home/bruce1996/nvme2/scRNA/GSE149614_coding_gene_hepatocyte_normalized.pkl")
#transform data.frame into seurat objectpbmc <- FindVariableFeatures(pbmc, selection.method = "vst")
pbmc <- CreateSeuratObject(counts = df, project = "GSE149614")
pbmc <- FindVariableFeatures(pbmc, selection.method = "vst",nfeatures=1346)
HVGs <- head(VariableFeatures(pbmc), 1346)

py_save_object(HVGs,filename = "/home/bruce1996/nvme2/scRNA/seurat_hvgs_list.pickle")

