import numpy as np
import pandas as pd
import pickle,random
from umap import UMAP
import hdbscan
import scanpy as sc
import anndata
from sklearn.metrics import normalized_mutual_info_score,adjusted_rand_score

def metadata_coverter(metadata) :
    converter = {}
    for idx,component in enumerate(sorted(metadata.unique())) :
        converter[component] = idx
    return converter

def calculate_umap_clustering_result(exp_matrix,ground_truth) :
    """_summary_

    Args:
        exp_matrix (pd.DataFrame): input expression matrix 
        ground_truth (np.array): metadata of samples
    """    
    reducer = UMAP(n_jobs=32,n_components=30)
    embedding = reducer.fit_transform(exp_matrix)
    hdbscan_labels = hdbscan.HDBSCAN().fit_predict(embedding)

    metric_dict = {}
    metric_dict['ARI'] = adjusted_rand_score(ground_truth,hdbscan_labels)
    metric_dict['NMI'] = normalized_mutual_info_score(ground_truth,hdbscan_labels)

    return metric_dict

def format_scanpy_input(exp_matrix,metadata) :
    """
    Args:
        exp_matrix (pd.DataFrame): gene expression matrix, row is sample ,col is gene 
        metadata (pd.DataFrame): metadata of sample in expression matrix
    """    
    x = exp_matrix.to_numpy()
    obs = metadata
    var = pd.DataFrame({'Gene' : exp_matrix.columns},index=exp_matrix.columns)
    data = anndata.AnnData(x,obs=obs,var=var)
    return data

def calculate_scanpy_leiden_clustering_result(exp_matrix,metadata,condition) :
    """
    Args:
        exp_matrix (pd.DataFrame): gene expression matrix, row is sample ,col is gene 
        metadata (pd.DataFrame): metadata of sample in expression matrix
        condition (str): colname of metadata
    """    
    data = format_scanpy_input(exp_matrix,metadata)
    stage_d = metadata_coverter(metadata[condition])
    ground_truth = metadata.loc[exp_matrix.index,condition].replace(stage_d).values
    
    sc.pp.neighbors(data, n_pcs=30,use_rep='X')
    sc.tl.umap(data)
    sc.tl.leiden(data, key_added="Leiden_res")

    metric_dict = {}
    metric_dict['ARI'] = adjusted_rand_score(ground_truth,data.obs.Leiden_res.values)
    metric_dict['NMI'] = normalized_mutual_info_score(ground_truth,data.obs.Leiden_res.values)

    return metric_dict

def main() :

    with open("/home/bruce1996/nvme2/scRNA/GSE149614_overall_normalized.pkl",'rb') as f :
        scRNA_exp_m = pickle.load(f)
    f.close()

    with open("/home/bruce1996/data/LIHC_anomaly_detection/data/coding_gene_info/ensembl2hgnc.pkl",'rb') as f :
        ens2hgsc = pickle.load(f)
    f.close()
    hgnc2ens = {}
    for ens,hgnc in ens2hgsc.items() :
        if isinstance(hgnc,str) :
            hgnc2ens[hgnc] = ens
    scRNA_metadata = pd.read_csv('/home/bruce1996/data/LIHC_anomaly_detection/validation_dataset/scRNA/GSE149614_HCC.metadata.txt',sep='\t',index_col=0)
    intersection = set(scRNA_exp_m.index).intersection(hgnc2ens.keys())
    ### subset expression matrix by coding gene & celltype
    exp_m = scRNA_exp_m.loc[intersection,scRNA_metadata.index[np.where((scRNA_metadata['Celltype'] == 'Hepatocyte') & (scRNA_metadata['site']=='Tumor'),True,False)]]
    exp_m.index = [hgnc2ens[x] for x in intersection]
    metadata = scRNA_metadata.loc[exp_m.columns,:]
    del scRNA_exp_m
    ### load vote information
    vote = pd.read_csv("/home/bruce1996/data/LIHC_anomaly_detection/manuscript/material/ensemble_learning_result/ensemble_hbv_only_np_ratio_35_vote_result.txt",sep='\t',index_col=0)
    voted_gene = list(set(vote.index[vote['Vote'] > 0]).intersection(exp_m.index))
    non_voted_gene = list(set(set(vote.index) - set(voted_gene)).intersection(exp_m.index))
    ### load hvg information
    seurat_hvg = pd.read_csv("/home/bruce1996/data/LIHC_anomaly_detection/manuscript/material/scRNA_seurat/seurat_hvg.txt",sep='\t',index_col=0)
    candidate_gene = list(set(seurat_hvg['EnsID'][np.where(seurat_hvg['EnsID'] != 'Non-coding',True,False)]).intersection(exp_m.index))
    non_candidate_gene = list(set(set(vote.index) - set(voted_gene)).intersection(exp_m.index))

    for condition in ['virus'] :
        stage_d = metadata_coverter(metadata[condition])
        y_true = metadata.loc[exp_m.columns,condition].replace(stage_d).values
        ### calculate random sample clustering result
        n_iteration = 100
        metric_df = pd.DataFrame(np.zeros([n_iteration,2]),columns=['ARI','NMI'])

        for idx in range(n_iteration) :
            print('%d st iteration for random sampling UMAP clustering' % (idx+1))
            candidate = random.sample(non_candidate_gene,len(candidate_gene))
            input = exp_m.loc[candidate,:]
            try :
                measurement = calculate_umap_clustering_result(input.T,y_true)
                metric_df.iloc[idx,0] = measurement['ARI'] 
                metric_df.iloc[idx,1] = measurement['NMI']
            except :
                pass
            '''
            leiden_measurement = calculate_scanpy_leiden_clustering_result(input.T,metadata,condition)
            leiden_df.iloc[idx,0] = leiden_measurement['ARI'] 
            leiden_df.iloc[idx,1] = leiden_measurement['NMI']
            '''
        ### calculate voted gene clustering result
        measurement = calculate_umap_clustering_result(exp_m.loc[voted_gene,:].T,y_true)
        metric_df = metric_df.append(measurement,ignore_index=True)
        '''
        leiden_measurement = calculate_scanpy_leiden_clustering_result(exp_m.loc[voted_gene,:].T,metadata,condition)
        leiden_df = leiden_df.append(leiden_measurement,ignore_index=True)
        '''

        metric_df.to_csv("/home/bruce1996/data/LIHC_anomaly_detection/manuscript/material/scRNA_measurement/%s_seurat_hvg_random_sampling_cluster_result.txt" % condition,sep='\t')
        #leiden_df.to_csv("/home/bruce1996/data/LIHC_anomaly_detection/manuscript/material/scRNA_measurement/%s_random_sampling_leiden_cluster_result.txt" % condition,sep='\t')

if __name__ == '__main__' :
    main()