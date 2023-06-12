library(monocle3)
library(ggplot2)
library(dplyr)
# load data
path="/home/bruce1996/data/LIHC_anomaly_detection/manifold_transformation/identify_marker_gene/"
fig_path="/home/bruce1996/data/LIHC_anomaly_detection/fig/scRNA_result/identify_marker_gene/"
condition='tumor_only'
exp_matrix <- read.table(paste0(path,"GSE149614_scRNA_vote_gene_exp_matrix_",condition,'.txt'),sep = '\t',header = T,row.names = 1)
cell_metadata <- read.table(paste0(path,"GSE149614_scRNA_metadata_",condition,'.txt'),sep = '\t',header = T,row.names = 1)
gene_annotation <- read.table(paste0(path,condition,"_gene_annotation.txt"),sep = '\t',header = T,row.names = 1)

cds <- new_cell_data_set(expression_data = as.matrix(exp_matrix),
                         cell_metadata = cell_metadata,
                         gene_metadata = gene_annotation)
# pre-processing data 
cds = preprocess_cds(cds,method = 'PCA',num_dim = 50,norm_method = 'none')
reduced_cds <- reduce_dimension(cds, reduction_method="UMAP",umap.fast_sgd=TRUE)
# identify marker gene 
marker_test_res <- top_markers(cds, group_cells_by="kmeans", 
                               reference_cells=10000, cores=32)

write.table(as.data.frame(marker_test_res),file = paste0(path,"GSE149614_scRNA_marker_gene_info_",condition,'.txt'),
            sep = '\t',quote = FALSE,row.names = FALSE)

top_specific_markers <- marker_test_res %>%
  filter(fraction_expressing >= 0.10) %>%
  group_by(cell_group) %>%
  top_n(1, pseudo_R2)

top_specific_marker_ids <- unique(top_specific_markers %>% pull(gene_id))

rowData(cds)$gene_short_name <- row.names(gene_annotation)
marker_plot = plot_genes_by_group(cds,
                    top_specific_marker_ids,
                    group_cells_by="kmeans",
                    max.size=5)

ggsave(paste0(fig_path,"marker_gene_info_",condition,'.png'),marker_plot,dpi = 300,width = 6,height = 8)
