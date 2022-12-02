import numpy as np
import pandas as pd
import pickle,json
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
        gene = list(set(go_genedict[go]).intersection(vote_result.index))
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