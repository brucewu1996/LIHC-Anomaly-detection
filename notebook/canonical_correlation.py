import numpy as np
import pandas as pd
import re,os
import random
import matplotlib.pyplot as plt
import seaborn as sns
from sklearn.preprocessing import StandardScaler
from sklearn.cross_decomposition import CCA

def voted_gene_list(gene_list,vote_number,threshold) :
    gene_l = np.array(gene_list)
    voted_idx =  vote_number > threshold
    return gene_l[voted_idx]

def cca_by_voted(non_hallmark_vote_result,hallmark_gene,exp_matrix,threshold) :

    voted_gene = voted_gene_list(non_hallmark_vote_result.index,non_hallmark_vote_result['Vote'].values,threshold)
    hallmark_m = exp_matrix.loc[hallmark_gene,:].T.to_numpy()
    voted_m = exp_matrix.loc[voted_gene,:].T.to_numpy()
    #calculated CCA
    n_comps = min(hallmark_m.shape[1], voted_m.shape[1])
    cca_score = np.zeros(2)
    cca = CCA(n_components=n_comps,scale=False,max_iter=2000)
    hallmark_t, voted_t = cca.fit_transform(hallmark_m,voted_m)
    cca_score[0] = np.corrcoef(hallmark_t[:,0],voted_t[:,0])[0,1]
    if n_comps > 1 :
        cca_score[1] = np.corrcoef(hallmark_t[:,1],voted_t[:,1])[0,1]
    return cca_score

def non_voted_random_cca(non_hallmark_vote_result,hallmark_gene,exp_matrix,threshold) :

    voted_gene = voted_gene_list(non_hallmark_vote_result.index,non_hallmark_vote_result['Vote'].values,threshold)
    hallmark_m = exp_matrix.loc[hallmark_gene,:].T.to_numpy()
    non_voted_gene = list(set(exp_matrix.index) - set(voted_gene) - set(hallmark_gene))
    target_non_voted_gene = random.sample(non_voted_gene,len(voted_gene))
    non_voted_m = exp_matrix.loc[target_non_voted_gene,:].T.to_numpy()
    #calculated CCA
    n_comps = min(hallmark_m.shape[1], non_voted_m.shape[1])
    cca_score = np.zeros(2)
    cca = CCA(n_components=n_comps,scale=False,max_iter=2000)
    hallmark_t, voted_t = cca.fit_transform(hallmark_m,non_voted_m)
    cca_score[0] = np.corrcoef(hallmark_t[:,0],voted_t[:,0])[0,1]
    if n_comps > 1 :
        cca_score[1] = np.corrcoef(hallmark_t[:,1],voted_t[:,1])[0,1]
    return cca_score

def main() :

    vote_result = pd.read_csv('/home/bruce1996/data/LIHC_anomaly_detection/ensemble_result/vote_result/normal_only_std_by_gene_vote_np_ratio_15.txt',sep='\t',index_col=0)
    coding_gene_idx = [bool(re.search('Synthetic',x)) == False for x in vote_result.index]# type: ignore    
    coding_gene = vote_result.index[coding_gene_idx]
    hallmark = pd.read_csv('/home/bruce1996/data/LIHC_anomaly_detection/data/hallmark_gene/hallmark_protein_coding_ens_id.txt',sep='\t')
    hallmark_gene = hallmark['EnsID'].values
    non_hallmark_gene = list(set(coding_gene) - set(hallmark_gene))

    exp_profile = pd.read_csv('/home/bruce1996/data/LIHC_anomaly_detection/data/exp_profile/LIHC_coding_gene_fpkm.txt',sep='\t',index_col=0)
    scaler = StandardScaler()
    std_x = scaler.fit_transform(exp_profile.T.to_numpy())
    exp_m = pd.DataFrame(std_x.T,index=exp_profile.index,columns=exp_profile.columns)
    ## CCA section
    threshold = np.arange(0,1000,50)
    n_iter_per_bin = 100
    cca_score = np.zeros([len(threshold) * n_iter_per_bin,2])
    n_vote_gene = np.zeros(len(threshold)* n_iter_per_bin)
    conditions = ['tumor_only','normal_only','hbv_only','nonhbv_only']
    np_ratio = [25,15,20,20]
    for c_idx,con in enumerate(conditions) :
        print("Random CCA for condition : %s " % con)
        vote_result = pd.read_csv('/home/bruce1996/data/LIHC_anomaly_detection/ensemble_result/vote_result/%s_std_by_gene_vote_np_ratio_%d.txt' % (con,np_ratio[c_idx]),sep='\t',index_col=0)
        non_hallmark_vote_result = vote_result.loc[non_hallmark_gene,:]
        for idx,th in enumerate(threshold) :
            voted_gene = voted_gene_list(non_hallmark_vote_result,non_hallmark_vote_result['Vote'].values,th)
            if len(voted_gene) == 0 :
                print("Without any nonhallmark gene vote number above %s" % th)
                break
            else :
                n_vote_gene[idx] = len(voted_gene)

            for n_iter in range(n_iter_per_bin) :
                cca_result = cca_by_voted(non_hallmark_vote_result,hallmark_gene,exp_m,th) 
                if len(cca_result) < 2 :
                    cca_score[n_iter_per_bin*idx+n_iter,0] = cca_result[0]
                else :
                    cca_score[n_iter_per_bin*idx+n_iter,:] = cca_result

            if idx % 5 == 0 and idx > 0 :
                print("CCA complete %d iterations !" % idx)

        threshold_array = np.zeros(len(threshold)* n_iter_per_bin)
        for idx,th in enumerate(threshold) :
            threshold_array[n_iter_per_bin*idx : n_iter_per_bin*(idx+1)] = th
        voted_gene_cca_result =pd.DataFrame({'n_voted_genes' : n_vote_gene,'CCA_coef1' : cca_score[:,0],'CCA_coef2' : cca_score[:,1],'Threshold' : threshold_array}) 
        voted_gene_cca_result.to_csv('/home/bruce1996/data/LIHC_anomaly_detection/voted_gene_cca/%s_random_cca_result.txt' % con,sep='\t')

if __name__ == '__main__' :
    main()