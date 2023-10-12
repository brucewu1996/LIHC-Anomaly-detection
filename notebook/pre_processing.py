import numpy as np
import pandas as pd
import pickle

def pre_processing(exp_m_path,sample_list,output_path) :

    #load expression matrix
    chunk = pd.read_csv(exp_m_path,sep='\t',index_col=0,chunksize=1000) 
    exp_m = pd.concat(chunk)
    target_m = exp_m.loc[:,sample_list]
    with open(output_path,'wb') as f :
        pickle.dump(target_m,f)
    del chunk,exp_m
    f.close()


def main() :
    metadata_path = '/home/bruce1996/data/LIHC_anomaly_detection/validation_dataset/scRNA/GSE149614_HCC.metadata.txt'
    exp_matrix_path = "/home/bruce1996/nvme2/scRNA/GSE149614_HCC.scRNAseq.S71915.normalized.txt"
    output_path = "/home/bruce1996/nvme2/scRNA/GSE149614_overall_normalized.pkl"

    metadata = pd.read_csv(metadata_path,sep='\t',index_col=0)
    sample_list = list(metadata.index)
    pre_processing(exp_matrix_path,sample_list,output_path)

if __name__ == '__main__' :
    main()