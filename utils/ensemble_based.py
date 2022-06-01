import pandas as pd
import numpy as np

def ensemble_result(result_path,threshold,ratio,n_run,n_sample,binary_flag = False) :

    vote = np.zeros(n_sample)

    for i in range(n_run) :
        file_name = result_path + 'lihc_ensemble_result_ratio_' + str(ratio) + '_test_' + str(i) + '.npy'
        result = np.load(file_name)
        vote = vote + result
    #transform vote number to binary 0/1
    if binary_flag :
            positive_idx = vote >= threshold * n_run
            vote[positive_idx] = 1
            vote[~positive_idx] = 0

    return vote

def gene_set_vote_number(gene_set,gene_list_of_ensemble_result,r_path,np_ratio,method = 'Mean') :
    
    vote_array = np.zeros([len(gene_set),len(np_ratio)])
    for idx,r in enumerate(np_ratio) :
        result = ensemble_result(r_path,0,r,10,19560)
        gene_idx = [x in gene_set for x in gene_list_of_ensemble_result]
        gene_vote = result[gene_idx]
        vote_array[:,idx] = gene_vote
    if method == 'Mean' :
        vr = np.mean(vote_array,axis=1)
    elif method == 'Median' :
        vr = np.median(vote_array,axis = 1)

    return vr   