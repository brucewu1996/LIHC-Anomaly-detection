import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import pickle
from scipy.stats import ranksums
from sklearn.cluster import KMeans
from statsmodels.sandbox.stats.multicomp import multipletests

class scRNA() :
    def __init__(self,exp_m,metadata,embedding,n_cluster) :
        self.expression_matrix = exp_m
        self.embedding = embedding
        self.genelist = list(exp_m.columns)
        self.sample = list(exp_m.index)
        self.metadata = metadata
        self.n_cluster = n_cluster
        self.foldchange = [np.zeros(len(self.genelist)) for _ in range(self.n_cluster) ] 
        #self.pvalue = [np.zeros(len(self.genelist)) for _ in range(self.n_cluster) ] 
        self.rejected = [np.zeros(len(self.genelist)) for _ in range(self.n_cluster) ] 
        self.adjusted_pvalue = [np.zeros(len(self.genelist)) for _ in range(self.n_cluster) ] 
        self.marker_gene_df = pd.DataFrame()
    '''
    def wilconox(self,condition,target,cluster_idx) :
        if condition not in self.metadata.columns :
            print("Make sure condition information in metadata")
            return
        candidate = self.metadata.index[np.where(self.metadata[condition] == target,True,False)]
        for idx,gene in enumerate(self.genelist) :
            _,pv = ranksums(self.expression_matrix.loc[candidate,gene].values,self.expression_matrix.loc[[x not in candidate for x in self.expression_matrix.index],gene].values)
            self.pvalue[cluster_idx][idx] = float(pv)
    '''
    def kmeans_clustering(self) :
        kmeans = KMeans(n_clusters=self.n_cluster, random_state=0).fit(self.embedding.loc[:,["UMAP1","UMAP2"]].to_numpy()) # type: ignore
        label = kmeans.labels_
        self.metadata['kmeans'] = label

    def bonferroni(self,condition,target,cluster_idx) :
        pv_list = []
        candidate = self.metadata.index[np.where(self.metadata[condition] == target,True,False)]
        for gene in self.genelist :
            _,pv = ranksums(self.expression_matrix.loc[candidate,gene].values,self.expression_matrix.loc[[x not in candidate for x in self.expression_matrix.index],gene].values)
            pv_list.append(pv)
        self.rejected[cluster_idx],self.adjusted_pvalue[cluster_idx], _, _ = multipletests(pv_list,alpha=0.05, method='bonferroni', is_sorted=False)

    def fold_change(self,condition,target,cluster_idx) :
        candidate = self.metadata.index[np.where(self.metadata[condition] == target,True,False)]
        for idx,gene in enumerate(self.genelist) :
            candidate_exp = self.expression_matrix.loc[candidate,gene].mean()
            no_candidate_exp = self.expression_matrix.loc[[x not in candidate for x in self.expression_matrix.index],gene].mean()
            if candidate_exp == 0  :
                self.foldchange[cluster_idx][idx] = np.inf
            elif no_candidate_exp == 0 :
                self.foldchange[cluster_idx][idx] = np.nan
            else :
                self.foldchange[cluster_idx][idx] = candidate_exp / no_candidate_exp
    def identify_marker_gene(self) :
        for i in range(self.n_cluster) :
            print("Identify marker gene for clsuter %d" % i)
            self.bonferroni('kmeans',i+1,i)
            self.fold_change('kmeans',i+1,i)
            
    def summarize_result(self) :
        m = np.empty([len(self.genelist),6],dtype=object)
        for i in range(self.n_cluster) :
            rejected = self.rejected[i]
            fold_change = self.foldchange[i]
            idx = fold_change > 2
            for g in range(len(self.genelist)) :
                if rejected[g] and idx[g] :
                    m[g,i] = 'Marker gene'
                else :
                    m[g,i] = 'Non marker gene'
        df = pd.DataFrame(m,index=self.genelist,columns=['Cluster' + str(x) for x in range(1,self.n_cluster+1)])
        self.marker_gene_df = df

def main() :
    metadata = pd.read_csv('/home/bruce1996/data/LIHC_anomaly_detection/validation_dataset/scRNA/GSE149614_HCC.metadata.txt',sep='\t',index_col=0)
    idx = np.where((metadata['site'] =='Tumor') & (metadata['Celltype'] =='Hepatocyte'),True,False)
    target = metadata.index[idx]
    target_metadata = metadata.loc[target,:]
    ##load expression information
    condition = ['normal_only','hbv_only','nonhbv_only']
    for con in condition :
            embedding_df = pd.read_csv("/home/bruce1996/data/LIHC_anomaly_detection/manifold_transformation/GSE149614_scRNA_umap_transformation_%s.txt" % con,sep='\t',index_col=0)
            sub_pickle_path = '/home/bruce1996/nvme2/scRNA/GSE149614_HCC_tumor_hepatocyte_%s.normalized.pickle' % con
            with open(sub_pickle_path,'rb') as f :
                    exp_m = pickle.load(f)
            con_scRNA = scRNA(exp_m=exp_m.astype(float),metadata=target_metadata,embedding=embedding_df,n_cluster=6)
            con_scRNA.kmeans_clustering()
            con_scRNA.identify_marker_gene()
            con_scRNA.summarize_result()
            df = con_scRNA.marker_gene_df
            df.to_csv("/home/bruce1996/data/LIHC_anomaly_detection/manifold_transformation/marker_gene_list_%s.txt" % con ,sep='\t')

if __name__ == '__main__' :
    main()