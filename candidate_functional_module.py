import numpy as np
import pandas as pd
import pickle,json
import traceback
import argparse
import matplotlib.pyplot as plt
import seaborn as sns

def candidate_go_vote_boxplot(go_genedict,vote_result,candidate_go,output_path,fig_title,fig_format) :
    '''
    go_genedict : dict ; key : go term ,value : list of gene belong to this GO term based on biomart annotation
    vote_result : dataframe ; ensemble learning vote result of each gene, index : gene, Columns "Vote" is vote number
    candidate_go : list ; list of candidate GO term
    '''
    plot_df = pd.DataFrame(columns=['GO','Gene','Vote'])
    vote_median = np.zeros(len(candidate_go.index))
    for go_idx,go in enumerate(candidate_go.index) :
        gene = list(set(go_genedict[go]).intersection(set(vote_result.index)))
        go_vote = vote_result.loc[gene,'Vote'].values
        vote_median[go_idx] = vote_result.loc[gene,'Vote'].median()
        df = pd.DataFrame({'GO' : [go] * len(gene),'Gene' : gene,'Vote' : go_vote})
        plot_df = pd.concat([plot_df,df])
        
    #import go description
    graph_d = dict()
    graph_path = '/home/bruce1996/data/GO/networkx/'
    go_field = ['biological_process','cellular_component','molecular_function']
    for field in go_field :
        file = graph_path + field + '.go'
        with open(file , 'rb') as f:  # notice the r instead of w
            graph_d[field] = pickle.load(f)
            f.close()

    with open("/home/bruce1996/data/GO/go2namespace.json",'rb') as f :
        go2namespace = json.load(f)
        f.close()
        
    order_df = pd.DataFrame({'GO' : list(candidate_go.index),'Vote' : vote_median}).sort_values(by='Vote',ascending=False)
    
    ### get go description
    go2description = dict()
    for go in candidate_go.index :
        namespace = go2namespace[go]
        go2description[go] = graph_d[namespace].nodes()[go]['name']
    plot_df['Description'] = plot_df['GO'].map(go2description)
    order_df['Description'] = order_df['GO'].map(go2description)
    n_go = candidate_go.shape[0]
    if n_go <= 20 : 
        plt.figure(figsize=(12,7))
    else :
        plt.figure(figsize=(round(n_go/3),7))
    sns.boxplot(data=plot_df,x='Description',y='Vote',palette='rainbow_r',order=list(order_df['Description']))
    plt.xticks(rotation=90)
    plt.xlabel('Candidate GO term')
    plt.title(fig_title)
    plt.savefig(output_path,dpi = 300,bbox_inches = 'tight',format = fig_format)
    
    return list(order_df['Description'])
    
def plot_module_performance(metric_df,output_path) :
    '''
    metric_df : dataframe; product by evaluate_functional_module_performance.py
    output_path : str; path of output fig
    '''
    
    plot_df = metric_df.drop(['Gene number'],axis=1).melt(id_vars='Description')
    plot_df.columns = ['Description','Metric','Value']
    if metric_df.shape[0] > 20 :
        plt.figure(figsize=(24,6))
    else :
        plt.figure(figsize=(12,6))
        
    sns.barplot(data=plot_df,x='Description',y='Value',hue='Metric',palette='Set2')
    plt.xticks(rotation=90)
    plt.xlabel("Functional module")
    plt.savefig(output_path,dpi=300,bbox_inches='tight')
    
def main() :
    parser = argparse.ArgumentParser()
    parser.add_argument("-g", "--go",help="path of candidate go term")
    parser.add_argument("-v", "--vote",help="path of vote result")
    parser.add_argument("-d", "--go_2_gene_dict",help="path of candidate go to gene dict")
    parser.add_argument("-l","--logistic",help = 'hbv prediction result by logistic regression')
    parser.add_argument("-r","--rf",help = 'hbv prediction result by random forest')
    parser.add_argument("-p","--prefix",help = 'prefix')
    parser.add_argument("-o","--output",help = 'fig output path')
    args = parser.parse_args()
    
    #This line opens a log file
    with open("log.txt", "a") as log:
        try:
            # some code
            print("plot candidate GO vote distribution", file = log)
            candidate_go = pd.read_csv(args.go,sep='\t',index_col=0)
            vote_result = pd.read_csv(args.vote,sep='\t',index_col=0)
            with open(args.go_2_gene_dict,'rb') as f :
                go_gene_dict = pickle.load(f)
            f.close()
            
            go_order = candidate_go_vote_boxplot(go_gene_dict,vote_result,candidate_go,args.output + 'Candidate_go_vote_' + args.prefix + '.png',"Candidate GO term vote distribution " + args.prefix,'png')
            logistic_result = pd.read_csv(args.logistic,sep='\t',index_col=0)
            randomforest_result = pd.read_csv(args.rf,sep='\t',index_col=0)
            
            gene_number = [int(x.split('=')[1][:-1]) for x in logistic_result['Description']]
            des = [x.split('(')[0][:-1] for x in logistic_result['Description']]
            des2genenumber = dict(zip(des,gene_number))
            plot_order = [x + " (n=%d)" % des2genenumber[x] for x in go_order]
            plot_module_performance(logistic_result,args.output + 'hbv_prediction_by_logistic_' + args.prefix + '.png')
            plot_module_performance(randomforest_result,args.output + 'hbv_prediction_by_randomforest_' + args.prefix + '.png')
        except Exception:
            traceback.print_exc(file=log)
            pass
    
if __name__ == '__main__' :
    main()