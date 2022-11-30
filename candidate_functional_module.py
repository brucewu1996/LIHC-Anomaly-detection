import numpy as np
import pandas as pd
import pickle
import os
import re
import matplotlib.pyplot as plt
import seaborn as sns

def candidate_go_vote_boxplot(go_genedict,vote_result,candidate_go,output_path,fig_title,fig_format) :
    '''
    go_genedict : dict ; key : go term ,value : list of gene belong to this GO term based on biomart annotation
    vote_result : dataframe ; ensemble learning vote result of each gene, index : gene, Columns "Vote" is vote number
    candidate_go : list ; list of candidate GO term
    '''
    plot_df = pd.DataFrame(columns=['GO','Gene','Vote'])
    for go in candidate_go.index :
        gene = list(set(go_genedict[go]).intersection(vote_result.index))
        go_vote = vote_result.loc[gene,'Vote'].values
        df = pd.DataFrame({'GO' : [go] * len(gene),'Gene' : gene,'Vote' : go_vote})
        plot_df = pd.concat([plot_df,df])
        
    plt.figure(figsize=(12,7))
    sns.boxplot(data=plot_df,x='GO',y='Vote',palette='rainbow_r')
    plt.xticks(rotation=90)
    plt.xlabel('Candidate GO term')
    plt.title(fig_title)
    plt.savefig(output_path,dpi = 300,bbox_inches = 'tight',format = fig_format)