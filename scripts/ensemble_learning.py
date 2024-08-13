import numpy as np
import pandas as pd
import random,os,pickle
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

class ensemble_learning_result:

    def __init__(self,hallmark_ens_list,hallmark_gene_symbol,gene_list):
        self.hallmark_gene = hallmark_ens_list
        self.hallmark_gene_symbol = hallmark_gene_symbol
        self.gene_list = gene_list
        self.hallmark_coversion_dict = dict(zip(self.hallmark_gene,self.hallmark_gene_symbol))
        self.svm_vote_result = None

    def generate_gene_vote_number(self) :
        if not self.svm_vote_result :
            return 
        self.vote_result = dict(zip(self.gene_list,self.svm_vote_result))


class svm_ensemble_learning() :
    def __init__(self,all_m,y_all,y_hallmark,r = 1,run = 1000,t = 20):
        self.overall_matrix = all_m
        self.y_overall = y_all
        self.y_hallmark = y_hallmark
        self.ratio = r
        self.run = run
        self.threads = t
        '''
        pos_idx = self.y_overall == 1
        hallmark_idx = self.y_hallmark == 1
        n_idx = self.y_overall == 0
        '''
        self.positive_matrix = self.overall_matrix[self.y_overall == 1,:]
        self.hallmark_matrix = self.overall_matrix[self.y_hallmark == 1,:]
        self.negative_matrix = self.overall_matrix[self.y_overall == 0,:]
        
    def single_svm(self,n) :
        population = range(self.negative_matrix.shape[0])
        # number of negative sample is required
        n_sample = self.ratio * self.positive_matrix.shape[0] 
        idx = random.sample(population,n_sample)
        negative = self.negative_matrix[idx,:]
        positive = self.positive_matrix
        X = np.concatenate((positive,negative))
        y = np.concatenate((np.repeat(1,positive.shape[0]),np.repeat(0,negative.shape[0]) ))

        x_train,x_test,y_train,y_test = train_test_split(X,y,test_size = 0.2,stratify = y)
        svm = SVC(kernel='linear')
        svm.fit(x_train,y_train)
        # svm predict gene is hallmark or not, it is vote number 
        y_pred = svm.predict(x_test)
        y_pred_overall = svm.predict(self.overall_matrix)
        vote = y_pred_overall
        # metric
        # precision = test dataset (20% of training)
        # precision_overall = entire expression matrix
        precision = precision_score(y_test,y_pred)
        precision_overall = precision_score(self.y_overall,y_pred_overall)
        recall = recall_score(y_test,y_pred)
        recall_overall = recall_score(self.y_overall,y_pred_overall)
        coef = svm.coef_[0]
        #hallmark metric 
        hallmark_x = np.concatenate((np.array(self.hallmark_matrix),np.array(negative)))
        hallmark_y = np.concatenate((np.repeat(1,self.hallmark_matrix.shape[0]),np.repeat(0,negative.shape[0])))
        hallmark_overall = np.concatenate((np.array(self.hallmark_matrix),np.array(self.negative_matrix)))
        hallmark_yoverall = np.concatenate((np.repeat(1,self.hallmark_matrix.shape[0]),np.repeat(0,self.negative_matrix.shape[0])))
        # hallmark = hallmark gene + training negative genes
        # hallmark overall = hallmark gene+ entire non-hallmark gene (not include the synthetic data)
        hallmark_pred_overall = svm.predict(hallmark_overall)
        hallmark_pred = svm.predict(hallmark_x)

        precision_hallmark = precision_score(hallmark_y,hallmark_pred)
        precision_overall_hallmark = precision_score(hallmark_yoverall,hallmark_pred_overall)
        recall_hallmark = recall_score(hallmark_y,hallmark_pred)
        recall_overall_hallmark = recall_score(hallmark_yoverall,hallmark_pred_overall)

        return vote,coef,precision,precision_overall,recall,recall_overall,precision_hallmark,precision_overall_hallmark,recall_hallmark,recall_overall_hallmark

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
    # set parameter 
    ROOT_DIR = '/home/bruce1996/LIHC-Anomaly-detection'
    EXP_M = f'{ROOT_DIR}/data/Ensemble-learning-training-data/lihc_coding_gene_std_by_gene_hbv_only_with_synthetic.txt'
    ENS_LIST = f'{ROOT_DIR}/data/Hallmark-information/hallmark_coding_gene.txt'
    STAND = False
    RUN = 50
    THREADS = 50
    OUTPUT_DIR = f'{ROOT_DIR}/data/Ensemble-learning-result/hbv_with_synthetic'
    PREFIX = 'hbv_only'
    ##
    hallmark_list = pd.read_csv(ENS_LIST,sep = '\t')
    exp_df = pd.read_csv(EXP_M,sep = '\t',index_col=0)
    if STAND :
        standard_exp = standardization(exp_df)
        standard_exp = standard_exp.T  # type: ignore
    else :
        standard_exp = exp_df.to_numpy()

    if os.path.isdir(OUTPUT_DIR) == False :
        os.mkdir(OUTPUT_DIR)
    ### subset by hallmark associated gene
    hallmark_ens_id = hallmark_list['Ensembl_ID'].values    
    hallmark_gene_symbol = hallmark_list['Gene_symbol'].values    
    hallmark_idx = [x in hallmark_ens_id for x in exp_df.index]
    syn_idx = [bool(search('Synthetic',x)) for x in exp_df.index]

    y_overall = np.zeros(exp_df.shape[0])
    y_overall[hallmark_idx] = 1
    y_overall[syn_idx] = 1
    y_hallmark = np.zeros(exp_df.shape[0])
    y_hallmark[hallmark_idx] = 1

    n_run = RUN
    n_thread = THREADS
    measurements = ['precision_test','precision_overall','recall_test','recall_overall',
                    'precision_hallmark','precision_hallmark_overall','recall_hallmark','recall_hallmark_overall']
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
            coef_r[i,:] = result[i,1]
        # a class record the ensemble learning output
        ensemble_record = ensemble_learning_result(hallmark_ens_list=hallmark_ens_id,
                                                hallmark_gene_symbol=hallmark_gene_symbol,
                                                gene_list=exp_df.index)

        ensemble_record.svm_vote_result = np.sum(vote_r,axis = 0)
        ensemble_record.coef = coef_r
        for m_idx,measurement in enumerate(measurements) :
            setattr(ensemble_record,measurement,result[:,m_idx+2])

        with open(f'{OUTPUT_DIR}/{PREFIX}_np_ratio_{r}_ensemble_result.pkl','wb') as f :
            pickle.dump(ensemble_record,f)
        f.close()
        del ensemble_record

if __name__ == '__main__':
    main()