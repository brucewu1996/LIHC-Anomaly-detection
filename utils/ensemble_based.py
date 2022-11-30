import pandas as pd
import numpy as np

def ensemble_result(result_path,prefix,threshold,ratio,binary_flag = False) :

    file_name = result_path + prefix + '_vote_np_ratio_'+str(ratio)+ '.npy'
    p_r = np.load(file_name)
    #transform vote number to binary 0/1
    if binary_flag :
        positive_idx = p_r >= threshold 
        p_r[positive_idx] = 1  # type: ignore
        p_r[~positive_idx] = 0  # type: ignore

    return p_r

def gene_set_vote_number(gene_set,gene_list_of_ensemble_result,r_path,r_prefix,np_ratio,method = 'Mean') :
    
    vote_array = np.zeros([len(gene_set),len(np_ratio)])
    for idx,r in enumerate(np_ratio) :
        result = ensemble_result(r_path,r_prefix,0,10)
        gene_idx = [x in gene_set for x in gene_list_of_ensemble_result]
        gene_vote = result[gene_idx]
        vote_array[:,idx] = gene_vote
    if method == 'Mean' :
        vr = np.mean(vote_array,axis=1)
    elif method == 'Median' :
        vr = np.median(vote_array,axis = 1)
    else :
        vr = vote_array

    return vr   

def merge_vote_result(r_path,r_prefix,index):
    
    vote_matrix = np.zeros([len(index),10])
    for r_idx,ratio in enumerate(range(5,55,5)) :
        vote = ensemble_result(r_path,r_prefix,0,ratio)
        vote_matrix[:,r_idx] = vote
    vote_df = pd.DataFrame(vote_matrix,index=index,columns=list(range(5,55,5)))
    
    return vote_df