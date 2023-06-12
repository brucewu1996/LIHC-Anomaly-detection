import pandas as pd
import numpy as np
import seaborn as sns
import matplotlib.pyplot as plt
from adjustText import adjust_text
from math import log10

def deseq2_volcano_plot(deseq2_result_path,fig_size = (10,7),fig_title='',fig_format='png',fig_output_path='test.png') :
    df = pd.read_csv(deseq2_result_path,sep='\t',index_col=0)
    df['logpvalue'] = df['pvalue'].apply(lambda x : 0 if x == 0 else -1 * log10(x))
    df['Class'] = np.where(((df['log2FoldChange'] >= 1)|(df['log2FoldChange'] <= -1)) & (df['pvalue'] < 0.05),"Significant","non-Significant")
    #volcano plot section
    plt.figure(figsize=fig_size)
    sns.scatterplot(data=df,x='log2FoldChange',y='logpvalue',hue='Class',palette='Set2')
    xlimit = max(abs(df['log2FoldChange'].min()),df['log2FoldChange'].max()) # type: ignore
    plt.xlim([xlimit+0.5,-1 * xlimit-0.5]) # type: ignore
    sig_genes = df.index[df['Class'] == 'Significant']
    #annotate differential expression genes
    xaxis = df.loc[sig_genes,'log2FoldChange'].values
    yaxis = df.loc[sig_genes,'logpvalue'].values
    texts = []
    for x, y, s in zip(xaxis,yaxis, sig_genes):
        texts.append(plt.text(x, y, s))
    adjust_text(texts, only_move={'points':'y', 'texts':'y'}, arrowprops=dict(arrowstyle="->", color='r', lw=0.5))
    #Add fig staff
    plt.ylabel("-log10pvalue")
    plt.title("%s (%d / %d)" % (fig_title,sum(df['Class'] == 'Significant'),df.shape[0]) )# type: ignore
    plt.savefig(fig_output_path,dpi=300,format=fig_format,bbox_inches='tight')

def deseq2_degenes(deseq2_result_path) :
    df = pd.read_csv(deseq2_result_path,sep='\t',index_col=0)
    df['Class'] = np.where(((df['log2FoldChange'] >= 1)|(df['log2FoldChange'] <= -1)) & (df['pvalue'] < 0.05),"Significant","non-Significant")
    sig_gene = list(df.index[df['Class'] == 'Significant'])
    return sig_gene