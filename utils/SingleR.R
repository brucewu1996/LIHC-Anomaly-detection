library(celldex)
library(scRNAseq)
library(SingleR)

matrix= read.table("/home/bruce1996/data/LIHC_anomaly_detection/validation_dataset/scRNA/GSE149614_HCC_hepatocyte_normalized.txt",header = T,row.names = 1,
                         sep = '\t')
ref.se <-celldex::HumanPrimaryCellAtlasData(ensembl=TRUE)
pred.hesc <- SingleR(test = m, ref = ref.se, assay.type.test='normalized',
                     labels = ref.se$label.fine)
cell_annotation = as.data.frame(pred.hesc$labels,row.names = colnames(m))
write.table(cell_annotation,file = "/home/bruce1996/data/LIHC_anomaly_detection/manifold_transformation/cell_annotation.txt",sep = '\t',quote = F)



