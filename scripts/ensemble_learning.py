import numpy as np
import pandas as pd
import random
import argparse
import os
from re import search
from sklearn.preprocessing import StandardScaler
from sklearn.model_selection import train_test_split
from sklearn.svm import SVC
from sklearn.metrics import confusion_matrix,precision_score,roc_auc_score,recall_score
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
    def __init__(self,all_m,y_all,y_hallmark,r = 1,run = 1000,t = 20):
        self.overall_matrix = all_m
        self.y_overall = y_all
        self.y_hallmark = y_hallmark
        self.ratio = r
        self.run = run
        self.threads = t
        pos_idx = self.y_overall == 1
        hallmark_idx = self.y_hallmark == 1
        n_idx = self.y_overall == 0
        self.positive_matrix = self.overall_matrix[pos_idx,:]
        self.hallmark_matrix = self.overall_matrix[hallmark_idx,:]
        self.negative_matrix = self.overall_matrix[n_idx,:]
        
    def single_svm(self,n) :
        population = range(self.negative_matrix.shape[0])
        n_sample = self.ratio * self.positive_matrix.shape[0]
        idx = random.sample(population,n_sample)
        negative = self.negative_matrix[idx,:]
        positive = self.positive_matrix
        X = np.concatenate((positive,negative))
        y = np.concatenate((np.repeat(1,positive.shape[0]),np.repeat(0,negative.shape[0]) ))

        x_train,x_test,y_train,y_test = train_test_split(X,y,test_size = 0.2,stratify = y)
        svm = SVC(kernel='linear')
        svm.fit(x_train,y_train)
        y_pred_overall = svm.predict(self.overall_matrix)
        y_pred = svm.predict(x_test)
        #metric
        precision = precision_score(y_test,y_pred)
        precision_overall = precision_score(self.y_overall,y_pred_overall)
        recall = recall_score(y_test,y_pred)
        recall_overall = recall_score(self.y_overall,y_pred_overall)
        coef = svm.coef_[0]
        #hallmark metric 
        hallmark = self.hallmark_matrix
        hallmark_n = self.negative_matrix
        hallmark_x = np.concatenate((np.array(hallmark),np.array(negative)))
        hallmark_y = np.concatenate((np.repeat(1,hallmark.shape[0]),np.repeat(0,negative.shape[0])))
        hallmark_overall = np.concatenate((np.array(hallmark),np.array(hallmark_n)))
        hallmark_yoverall = np.concatenate((np.repeat(1,hallmark.shape[0]),np.repeat(0,hallmark_n.shape[0])))

        hallmark_pred_overall = svm.predict(hallmark_overall)
        hallmark_pred = svm.predict(hallmark_x)

        precision_hallmark = precision_score(hallmark_y,hallmark_pred)
        precision_overall_hallmark = precision_score(hallmark_yoverall,hallmark_pred_overall)
        recall_hallmark = recall_score(hallmark_y,hallmark_pred)
        recall_overall_hallmark = recall_score(hallmark_yoverall,hallmark_pred_overall)

        return y_pred_overall,precision,precision_overall,recall,recall_overall,coef,precision_hallmark,precision_overall_hallmark,recall_hallmark,recall_overall_hallmark

    def ensemble_svm(self) :
        pool = mp.Pool(self.threads)
        n_run = self.run
        result = pool.map(self.single_svm,range(n_run))
        pool.close()
        pool.join()
        #merge parallel result
        metric_array = np.stack(result,axis = 0) # type: ignore
        
        return metric_array

def vote_for_ensemble(result,threshold) :
    vote = np.sum(result,axis = 0)
    vote[vote >= threshold] = 1  # type: ignore
    vote[vote < threshold] = 0  # type: ignore
    return vote

def main() :

    parser = argparse.ArgumentParser()
    parser.add_argument("-e", "--input",help="path of expression profile")
    #parser.add_argument("-s", "--synthetic",help="path of synthetic data")
    parser.add_argument("-c", "--ens_list",help="path of hallmark gene ens id")
    parser.add_argument("-r", "--run",type=int,help="number of svm for an ensemble model")
    parser.add_argument("-t", "--threads",type=int,help="number of threads to use")
    parser.add_argument("-o","--output_path",help = 'path of ensemble model output')
    parser.add_argument("-p","--prefix",type=str,help="prefix of metric output file")
    parser.add_argument("--stand",type=bool,default=False,help = "Standardziation of not,if standardzied choose False")
    args = parser.parse_args()

    hallmark_list = pd.read_csv(args.ens_list,sep = '\t')
    exp_df = pd.read_csv(args.input,sep = '\t',index_col=0)
    #syn_df = pd.read_csv(args.synthetic,sep = '\t',index_col=0)
    #exp_df = pd.concat([exp_df,syn_df],axis=0)
    if args.stand :
        standard_exp = standardization(exp_df)
        standard_exp = standard_exp.T  # type: ignore
    else :
        standard_exp = exp_df.to_numpy()

    if os.path.isdir(args.output_path) == False :
        os.mkdir(args.output_path)
    ### subset by hallmark associated gene
    ens_id = hallmark_list['EnsID'].values       
    hallmark_idx = [x in ens_id for x in exp_df.index]
    syn_idx = [bool(search('Synthetic',x)) for x in exp_df.index]

    y_overall = np.zeros(exp_df.shape[0])
    y_overall[hallmark_idx] = 1
    y_overall[syn_idx] = 1
    y_hallmark = np.zeros(exp_df.shape[0])
    y_hallmark[hallmark_idx] = 1

    n_run = args.run
    n_thread = args.threads
    ensemble = svm_ensemble_learning(standard_exp,y_overall,y_hallmark,run = n_run,t = n_thread)
    ###differ positive / negative ratio
    for r in np.arange(5,55,5) :
        print("Ensemble learning for N/P ratio  : %d " % r)
        ensemble.ratio = r
        result = ensemble.ensemble_svm()
        # result
        vote_r = np.zeros((n_run,exp_df.shape[0]))
        coef_r = np.zeros((n_run,exp_df.shape[1]))
        for i in range(n_run) :
            vote_r[i,:] = result[i,0]  # type: ignore
            coef_r[i,:] = result[i,5]
        
        precision_r = result[:,1]
        precision_overall_r = result[:,2]
        recall_r = result[:,3]
        recall_overall_r = result[:,4]

        precision_h = result[:,6]
        precision_overall_h = result[:,7]
        recall_h = result[:,8]
        recall_overall_h = result[:,9]

        vote_result = np.sum(vote_r,axis = 0)
        np.save(args.output_path + args.prefix +'_vote_np_ratio_'+str(r)+ '.npy',vote_result)
        np.save(args.output_path + args.prefix +'_precision_np_ratio_'+str(r)+ '.npy',precision_r)
        np.save(args.output_path + args.prefix +'_precision_overall_np_ratio_'+str(r)+ '.npy',precision_overall_r)
        np.save(args.output_path + args.prefix + '_recall_np_ratio_'+str(r)+ '.npy',recall_r)
        np.save(args.output_path + args.prefix + '_recall_overall_np_ratio_'+str(r)+ '.npy',recall_overall_r)
        np.save(args.output_path + args.prefix + '_coef_np_ratio_'+str(r)+ '.npy',coef_r)
        #hallmark
        np.save(args.output_path + args.prefix +'_precision_hallmark_np_ratio_'+str(r)+ '.npy',precision_r)
        np.save(args.output_path + args.prefix +'_precision_overall_hallmark_np_ratio_'+str(r)+ '.npy',precision_overall_r)
        np.save(args.output_path + args.prefix + '_recall_hallmark_np_ratio_'+str(r)+ '.npy',recall_r)
        np.save(args.output_path + args.prefix + '_recall_overall_hallmark_np_ratio_'+str(r)+ '.npy',recall_overall_r)

if __name__ == '__main__':
    main()