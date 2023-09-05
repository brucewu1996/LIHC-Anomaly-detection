import numpy as np
import pandas as pd
import pickle,random
from umap import UMAP
from sklearn.metrics import normalized_mutual_info_score,adjusted_rand_score,silhouette_score
from sklearn.cluster import KMeans

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
    n_cluster = len(set(ground_truth))
    reducer = UMAP(n_jobs=32)
    embedding = reducer.fit_transform(exp_matrix)
    kmeans_result = KMeans(n_clusters=n_cluster).fit_predict(embedding) # type: ignore

    metric_dict = {}
    metric_dict['ARI'] = adjusted_rand_score(ground_truth,kmeans_result)
    metric_dict['NMI'] = normalized_mutual_info_score(ground_truth,kmeans_result)
    metric_dict['SC'] = silhouette_score(embedding,labels=kmeans_result) # type: ignore

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

    condition = 'patient'
    stage_d = metadata_coverter(metadata[condition])
    y_true = metadata.loc[exp_m.columns,condition].replace(stage_d).values
    ### calculate random sample clustering result
    n_iteration = 100
    metric_df = pd.DataFrame(np.zeros([n_iteration,3]),columns=['ARI','NMI','SC'])

    for idx in range(n_iteration) :
        print('%d st iteration for random sampling UMAP clustering' % idx)
        candidate = random.sample(non_voted_gene,len(voted_gene))
        input = exp_m.loc[candidate,:]
        measurement = calculate_umap_clustering_result(input.T,y_true)
        metric_df.iloc[idx,0] = measurement['ARI'] 
        metric_df.iloc[idx,1] = measurement['NMI']
        metric_df.iloc[idx,2] = measurement['SC']

    metric_df.to_csv("/home/bruce1996/data/LIHC_anomaly_detection/manuscript/material/scRNA_measurement/%s_random_sampling_cluster_result.txt" % condition,sep='\t')

if __name__ == '__main__' :
    main()