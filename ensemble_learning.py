import numpy as np
import pandas as pd
import random
import argparse
from re import search
from sklearn.preprocessing import StandardScaler
from sklearn.model_selection import train_test_split
from sklearn.svm import SVC
from sklearn.metrics import confusion_matrix,precision_score
import multiprocessing as mp


def standardization(df) :
    '''
    df : exp profile {gene,sample}
    '''
    array = df.T.to_numpy()
    scaler = StandardScaler().fit(array)
    standardized_array = scaler.transform(array)

    return standardized_array

def npv(y_true,y_pred) :
    tn, _, fn, _ = confusion_matrix(y_true, y_pred).ravel()
    npv = tn / (tn+fn)
    return npv

class svm_ensemble_learning() :
    def __init__(self,all_m,y_all,r = 1,run = 1000,t = 20):
        self.overall_matrix = all_m
        self.y_overall = y_all
        self.ratio = r
        self.run = run
        self.threads = t
        pos_idx = self.y_overall == 1
        n_idx = self.y_overall == 0
        self.positive_matrix = self.overall_matrix[pos_idx,:]
        self.negative_matrix = self.overall_matrix[n_idx,:]

    def data_split(self) :
        '''
        exp_df : data frame of expression profile
        y : numpy array ,annotation of exp_df, must be binary
        '''
        exp_df = self.overall_matrix
        y = self.y_overall
        pos_idx = y == 1
        n_idx = y == 0
        self.positive_matrix = exp_df[pos_idx,:]
        self.negative_matrix = exp_df[n_idx,:]
        
        
    def single_svm(self,n) :
        population = range(self.negative_matrix.shape[0])
        n_sample = self.ratio * self.positive_matrix.shape[0]
        idx = random.sample(population,n_sample)
        negative = self.negative_matrix[idx,:]
        positive = self.positive_matrix
        X = np.concatenate((positive,negative))
        y = np.concatenate((np.repeat(1,positive.shape[0]),np.repeat(0,negative.shape[0]) ))

        x_train,x_test,y_train,y_test = train_test_split(X,y,test_size = 0.2)
        svm = SVC(kernel='linear')
        svm.fit(x_train,y_train)
        y_pred_overall = svm.predict(self.overall_matrix)
        y_pred = svm.predict(x_test)
        precision = precision_score(y_test,y_pred)

        return y_pred_overall,precision

    def ensemble_svm(self) :
        pool = mp.Pool(self.threads)
        n_run = self.run
        result = pool.map(self.single_svm,range(n_run))
        pool.close()
        pool.join()
        #merge parallel result
        metric_array = np.stack(result,axis = 0)
        
        return metric_array

def vote_for_ensemble(result,threshold) :
    vote = np.sum(result,axis = 0)
    vote[vote >= threshold] = 1
    vote[vote < threshold] = 0
    return vote

def main() :

    parser = argparse.ArgumentParser()
    parser.add_argument("-i", "--input",help="path of lihc expression profile")
    parser.add_argument("-c", "--conversion_table",help="path of hbv conversion table")
    parser.add_argument("-r", "--run",type=int,help="number of svm for an ensemble model")
    parser.add_argument("-t", "--threads",type=int,help="number of threads to use")
    parser.add_argument("-o","--output_path",help = 'path of ensemble model output')
    parser.add_argument("-p","--prefix",type=str,help="prefix of metric output file")
    parser.add_argument("-s","--stand",type=bool,default=False,help = "Standardziation of not,if standardzied choose False")
    args = parser.parse_args()

    hbv_conversion_table = pd.read_csv(args.conversion_table,sep = '\t',index_col=0)
    exp_df = pd.read_csv(args.input,sep = '\t',index_col=0)
    if args.stand :
        standard_exp = standardization(exp_df)
        standard_exp = standard_exp.T
    else :
        standard_exp = exp_df.to_numpy()
    ### subset by hbv associated gene
    ens_id = hbv_conversion_table['Ensembl_ID'].values       
    hbv_idx = [x.split('.')[0] in ens_id for x in exp_df.index]
    syn_idx = [bool(search('Synthesized',x)) for x in exp_df.index]

    y_overall = np.zeros(exp_df.shape[0])
    y_overall[hbv_idx] = 1
    y_overall[syn_idx] = 1

    n_run = args.run
    n_thread = args.threads
    ensemble = svm_ensemble_learning(standard_exp,y_overall,run = n_run,t = n_thread)
    ###differ positive / negative ratio
    for r in np.arange(5,55,5) :
        ensemble.ratio = r
        result = ensemble.ensemble_svm()
        # result
        vote_r = np.zeros((n_run,exp_df.shape[0]))
        precision_r = np.zeros(n_run)
        for i in range(n_run) :
            vote_r[i,:] = result[i,0]
            precision_r[i] = result[i,1]
        vote_result = np.sum(vote_r,axis = 0)
        np.save(args.output_path + args.prefix +'_vote_np_ratio_'+str(r)+ '.npy',vote_result)
        np.save(args.output_path + args.prefix +'_precision_np_ratio_'+str(r)+ '.npy',precision_r)

if __name__ == '__main__':
    main()