import datetime
import numpy as np
import pandas as pd
import os,sys,re
import matplotlib.pyplot as plt
import seaborn as sns
sys.path.append("/home/bruce1996/repo/LIHC_anomaly_detection/")
from sklearn.model_selection import cross_val_score
from sklearn.svm import SVC
from joblib import parallel_backend
from ensemble_model_performance import *

parallel_backend(backend='threading', n_jobs=12)

def read_vote_result(result_path,prefix,np_ratio,idx) :
    vote = np.load(result_path + '%s/%s__vote_np_ratio_%d.npy' % (prefix,prefix,np_ratio))
    df = pd.Series(vote,index=idx)
    idx = df.index.str.contains("ENSG")
    return df[idx]

def read_ensemble_result(result_path,prefix,np_ratio,gene_idx) : 
    """
    Args:
        result_path (str): path of ensemble learning model result
        prefix (str): prefix of ensemble learning model
        np_ratio (int): N/P ratio of ensemble learning model
        gene_idx (list): Gene index of ensemble learning model input

    Returns:
        df[idx]: return metric of model 
        vote[idx] : return vote number without synthetic data
    """    

    vote = np.load(result_path + '%s/%s_vote_np_ratio_%d.npy' % (prefix,prefix,np_ratio))
    precision = np.load(result_path + '%s/%s_precision_np_ratio_%d.npy' % (prefix,prefix,np_ratio),allow_pickle=True)
    overall_p = np.load(result_path + '%s/%s_precision_overall_np_ratio_%d.npy' % (prefix,prefix,np_ratio),allow_pickle=True)
    recall = np.load(result_path + '%s/%s_recall_np_ratio_%d.npy' % (prefix,prefix,np_ratio),allow_pickle=True)
    overall_r = np.load(result_path + '%s/%s_recall_overall_np_ratio_%d.npy' % (prefix,prefix,np_ratio),allow_pickle=True)
    df = pd.DataFrame({'Precision' : precision,'Precision_overall' : overall_p,
                       'Recall' : recall, 'Recall_overall' : overall_r})
    vote = pd.Series(vote,index=gene_idx)
    idx = vote.index.str.contains("ENSG")
    return df,vote[idx]

def model_performance(result_path,prefix,hallmark,gene_idx,metric,output_path,format='png') :

    no_syn_idx = gene_idx.str.contains("ENSG")
    metric_median = np.zeros(10)
    metric_std = np.zeros(10)
    overall_metric_median = np.zeros(10)
    overall_metric_std = np.zeros(10)
    auc_array = np.zeros(10) 
    auc_std_array = np.zeros(10)
    for np_idx,np_r in enumerate(np.arange(5,55,5)) :
        start = datetime.datetime.now()
        print("%s Predict hallmark gene by vote number of np ratio %d in %s condition" % (start.strftime("%H:%M:%S"),np_r,prefix))  
        metric_df,vote_df = read_ensemble_result(result_path,"%s_std_by_gene" % prefix,np_r,gene_idx)
        metric_median[np_idx] = np.median(metric_df[metric].values)
        metric_std[np_idx] = np.std(metric_df[metric].values)
        overall_metric_median[np_idx] = np.median(metric_df[metric+'_overall'].values)
        overall_metric_std[np_idx] = np.std(metric_df[metric+'_overall'].values)
        no_syn_gene = gene_idx[no_syn_idx]
        x = np.array(vote_df[no_syn_gene].values)
        y = np.array([1 if x in hallmark else 0 for x in no_syn_gene])
        
        #x_train,x_test,y_train,y_test = train_test_split(x.reshape(-1, 1),y.reshape(-1, 1),test_size = 0.2,stratify = y)
        svm = SVC(kernel='linear')
        cv_score = cross_val_score(svm,x.reshape(-1,1),y,cv=5,scoring='roc_auc')
        auc_array[np_idx] = np.median(cv_score)
        auc_std_array[np_idx] = cv_score.std()
        end = datetime.datetime.now()
        gap = end - start
        print("%s Predict hallmark gene by vote number of np ratio %d in %s condition completed, executime time : %s" % (end.strftime("%H:%M:%S"),np_r,prefix,str(gap).split('.')[0] + ' seconds')) 
        #fpr, tpr, threshold = roc_curve(y, x, pos_label=1)
        #auc_array[np_idx] = auc(fpr, tpr)

    plt.figure(figsize=(10,5))
    plt.errorbar(np.arange(5,55,5),metric_median,yerr = metric_std,marker = '*', label=metric.capitalize(),color = "darksalmon",ms = 10)
    plt.errorbar(np.arange(5,55,5),overall_metric_median,yerr = overall_metric_std, marker = '*',label='Overall ' +metric ,color = "mediumaquamarine",ms = 10)
    plt.errorbar(np.arange(5,55,5),auc_array,yerr=auc_std_array,marker = '*', label='AUC',color = "#FFAACF",ms = 10)
    plt.ylabel("Metrics")
    plt.ylim([0,1])
    plt.xlabel('N/P ratio')
    plt.legend()
    plt.title("Ensemble learning model performance (%s)" % (prefix))
    plt.savefig(output_path + 'ensemble_performance_predict_by_vote_%s_%s.%s' % (metric,prefix,format),bbox_inches = 'tight',dpi=300,format=format)

def main() :
    vote_path = '/home/bruce1996/data/LIHC_anomaly_detection/ensemble_result/without_synthetic/'
    training_path = '/home/bruce1996/data/LIHC_anomaly_detection/ensemble_training/'
    hallmark = pd.read_csv('/home/bruce1996/data/LIHC_anomaly_detection/data/hallmark_gene/hallmark_protein_coding_ens_id.txt',sep='\t')
    hallmark_genes = hallmark['EnsID'].values
    output_path = "/home/bruce1996/data/LIHC_anomaly_detection/manuscript/Fig1_identify_hallmark_gene/"
    ##
    prefix_list = ['hbv_only','tumor_only']
    for prefix in prefix_list :
        exp_m = pd.read_csv(training_path + "lihc_coding_gene_std_by_gene_%s.txt" % prefix,sep='\t',index_col=0)
        model_performance(vote_path,prefix,hallmark_genes,exp_m.index,'Precision',output_path,format='svg')


if __name__ == '__main__' :
    main()