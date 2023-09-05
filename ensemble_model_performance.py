import numpy as np
import argparse,os
import matplotlib.pyplot as plt
from utils.ensemble_based import ensemble_result

def metric_performance(path,prefix,metric,title,output_path,hallmark=False) :
    '''
    path : str; path of ensemble result
    prefix : str; prefix of ensemble result
    metric : str; metric to evaluate of ensemble result performance 
    title : str; figure title
    output_path : str; figure output path
    hallmark : boolean; evaluation ensemble result without any synthetic data
    '''
    
    metric_median = np.zeros(10)
    metric_std = np.zeros(10)
    for idx,r in enumerate(np.arange(5,55,5)) :
        if hallmark :
            tet = np.load(path + prefix + '_' + metric + '_hallmark_np_ratio_'+str(r)+ '.npy',allow_pickle=True)
        else :
            tet = np.load(path + prefix + '_' + metric + '_np_ratio_'+str(r)+ '.npy',allow_pickle=True)
        metric_median[idx] = round(np.median(tet),3)
        metric_std[idx] = round(np.std(tet),3)

    overall_metric_median = np.zeros(10)
    overall_metric_std = np.zeros(10)
    for idx,r in enumerate(np.arange(5,55,5)) :
        if hallmark :
            tet = np.load(path + prefix + '_' + metric + '_overall_hallmark_np_ratio_'+str(r)+ '.npy',allow_pickle=True)
        else :
            tet = np.load(path + prefix + '_' + metric + '_overall_np_ratio_'+str(r)+ '.npy',allow_pickle=True)
        overall_metric_median[idx] = round(np.median(tet),3)
        overall_metric_std[idx] = round(np.std(tet),3)


    plt.figure(figsize=(10,5))
    plt.errorbar(np.arange(5,55,5),metric_median,yerr = metric_std,marker = '*', label=metric.capitalize(),color = "darksalmon",ms = 10)
    plt.errorbar(np.arange(5,55,5),overall_metric_median,yerr = overall_metric_std, marker = '*',label='Overall ' +metric ,color = "mediumaquamarine",ms = 10)
    plt.ylabel(metric.capitalize())
    plt.xlabel('N/P ratio')
    plt.legend()
    plt.title(title)
    plt.savefig(output_path,bbox_inches = 'tight',dpi=300)

def main() :
    parser = argparse.ArgumentParser()
    parser.add_argument("-v", "--vote",help="path of vote result")
    parser.add_argument("-p", "--prefix",help="prefix of vote result")
    parser.add_argument("-o","--output_path",type=str,help = 'path of ensemble model output')
    parser.add_argument("-f","--figure_prefix",help = 'prefix of figure')
    parser.add_argument("-t","--title",type=str,default="Median performance of single SVM in ensemble learning model" )
    args = parser.parse_args()  

    ### plot
    result_path = args.vote
    result_prefix = args.prefix
    title = args.title
    
    if os.path.exists(args.output_path) == False :
        os.mkdir(args.output_path)
        
    metric_performance(result_path ,result_prefix,'recall',title,args.output_path + args.figure_prefix + 'ensemble_model_median_recall.png',hallmark=True)
    metric_performance(result_path ,result_prefix,'precision',title,args.output_path + args.figure_prefix + 'ensemble_model_median_precision.png',hallmark=True)
    
if __name__ == '__main__' :
    main()