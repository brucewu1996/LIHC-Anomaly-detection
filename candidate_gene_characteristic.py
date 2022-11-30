import numpy as np
import pandas as pd
import re,os,argparse
import matplotlib.pyplot as plt
from sklearn.metrics import recall_score

def ensemble_result(result_path,prefix,threshold,ratio,binary_flag = False) :

    file_name = result_path + prefix + '_vote_np_ratio_'+str(ratio)+ '.npy'
    p_r = np.load(file_name)
    #transform vote number to binary 0/1
    if binary_flag :
        positive_idx = p_r >= threshold 
        p_r[positive_idx] = 1  # type: ignore
        p_r[~positive_idx] = 0  # type: ignore

    return p_r


def candidate_hallmark_recall(hallmark_gene,exp_profile,result_path,result_prefix,np_ratio,fig_output,vote_floor=1000) :
    '''
    hallmark gene : list; Record hallmark gene EnsID
    exp_profile : dataframe; gene expression profile after standardzation , row is gene column is sample
    result path : str; path of ensemble learning result
    result prefix : str; prefix of ensemble learning result
    np_ratio : int; Value of ensemble learning NP ratio
    fig_output : str; path of fig output
    '''
    
    threshold = np.arange(1,vote_floor+1)
    pos_num = np.zeros(len(threshold))
    ens_vote_matrix = np.zeros((len(hallmark_gene),len(threshold)))
    ens_idx = [x in hallmark_gene for x in exp_profile.index]
    no_synthetic_idx = [bool(re.search('Synthetic',x)) == False for x in exp_profile.index]

    for idx,t in enumerate(threshold) :
        y_pred = ensemble_result(result_path,result_prefix,t,np_ratio,binary_flag=True)
        ens_vote_matrix[:,idx] = y_pred[ens_idx] 
        pos_num[idx] = sum(y_pred[no_synthetic_idx])  # type: ignore

    ens_pos_num = ens_vote_matrix.sum(axis=0)
    ens_ratio = ens_pos_num / pos_num

    ens_recall = np.zeros(len(threshold))
    ens_y = np.repeat(1,len(hallmark_gene))
    for i in range(len(threshold)) :
        ens_pred = ens_vote_matrix[:,i]
        ens_recall[i] = recall_score(ens_y,ens_pred)
        
    ###plot
    fig,ax1 = plt.subplots(figsize = (8,5))
    ax2 = plt.twinx(ax1)
    c1 = ax1.errorbar(threshold,ens_ratio,label = 'Hallmark ratio in candidate gene',color = 'darksalmon')
    c2 = ax1.errorbar(threshold,ens_recall,label = 'Hallmark recall',color = "aquamarine")
    ax1.set_ylabel("Metric")
    c3 = ax2.errorbar(threshold,pos_num,label = 'Number of positive',color = '#33CCFF')
    ax2.set_ylabel("Number of positive")

    curves = [c1,c2,c3]
    ax1.legend(curves,['Hallmark ratio in candidate gene','Hallmark recall','Number of positive'],loc = 'center left')
    plt.xlabel("Vote threshold")
    plt.title("Hallmark gene racall & ratio in candidate gene")
    plt.savefig(fig_output,bbox_inches = 'tight',dpi = 300)
    
def candidate_gene_proportion(hallmark_gene,exp_profile,result_path,result_prefix,np_ratio,fig_output,vote_floor=1000) :
    '''
    hallmark gene : list; Record hallmark gene EnsID
    exp_profile : dataframe; gene expression profile after standardzation , row is gene column is sample
    result path : str; path of ensemble learning result
    result prefix : str; prefix of ensemble learning result
    np_ratio : int; Value of ensemble learning NP ratio
    fig_output : str; path of fig output
    '''
    threshold = np.arange(1,vote_floor+1)
    pos_num = np.zeros(len(threshold))
    ens_pos = np.zeros(len(threshold))
    syn_pos = np.zeros(len(threshold))
    ens_idx = [x in hallmark_gene for x in exp_profile.index]
    synthetic_idx = [bool(re.search('Synthetic',x)) for x in exp_profile.index]

    for idx,t in enumerate(threshold) :
        y_pred = ensemble_result(result_path,result_prefix,t,np_ratio,binary_flag=True)
        pos_num[idx] = sum(y_pred)  # type: ignore
        ens_pos[idx] = sum(y_pred[ens_idx])  # type: ignore
        syn_pos[idx] = sum(y_pred[synthetic_idx])  # type: ignore

    residual_pos = (pos_num - ens_pos - syn_pos) / pos_num
    ens_pos = ens_pos/ pos_num
    syn_pos = syn_pos / pos_num

    df = pd.DataFrame({'Hallmark positive' : ens_pos,'Synthetic positive' : syn_pos,'Residual' : residual_pos})
    
    df.columns = ['Hallmark','Synthetic','Candidate']  # type: ignore
    df.plot(kind = 'bar',stacked = True,figsize = (12,5))
    plt.xticks(np.arange(0,vote_floor,100))
    plt.xlabel("Vote threshold",fontsize = 12)
    plt.ylabel("Ratio",fontsize =12)
    plt.title("Gene vote above threhold proportion")
    plt.legend(loc = 'center right')
    plt.savefig(fig_output,bbox_inches = 'tight',dpi = 300)
    
def main() :
    parser = argparse.ArgumentParser()
    parser.add_argument("-i", "--input",help="path of  expression profile")
    parser.add_argument("-r", "--reference",help="path of hallmark gene list")
    parser.add_argument("-v", "--vote",help="path of vote result")
    parser.add_argument("-p", "--prefix",help="prefix of vote result")
    parser.add_argument("-n", "--np_ratio",type=int,help="NP ratio")
    parser.add_argument("-o","--output_path",type=str,help = 'path of ensemble model output')
    parser.add_argument("-f","--figure_prefix",help = 'prefix of figure')
    args = parser.parse_args()
    
    exp_profile = pd.read_csv(args.input,sep='\t',index_col=0)
    
    hallmark = pd.read_csv(args.reference,sep='\t')
    hallmark.columns = ["EnsID"]  # type: ignore
    hallmark_gene = hallmark['EnsID'].values
    
    result_path = args.vote
    result_prefix = args.prefix
    if os.path.exists(args.output_path) == False:
        os.mkdir(args.output_path)   
        
    candidate_hallmark_recall(hallmark_gene,exp_profile,result_path,result_prefix,args.np_ratio,args.output_path +args.figure_prefix +'Candidate_hallmark_recall.png')
    candidate_gene_proportion(hallmark_gene,exp_profile,result_path,result_prefix,args.np_ratio,args.output_path +args.figure_prefix +'Candidate_composition.png')   

if __name__ == '__main__':
    main()