import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from sklearn.decomposition import PCA
import seaborn as sns
import argparse


def plot_pc_explain_ratio(exp_matrix,output_path,explained_ratio_threshold=0.8) :
    '''
    exp_matrix : numpy array; row is gene column is sample
    explained_ratio_threshold : float; threshold of sum of principle component explained ratio 
    '''
    sum_variance = [] 
    for i in range(2,100) :
        pca = PCA(n_components=i)
        exp_pca = pca.fit_transform(exp_matrix)
        explained_ratio = sum(pca.explained_variance_ratio_)
        sum_variance.append(explained_ratio)
        if explained_ratio > explained_ratio_threshold :
            break
        
    df = pd.DataFrame({'Explained variance' : sum_variance + list(np.repeat(explained_ratio_threshold,len(sum_variance))),
                       'X' : np.concatenate([np.arange(len(sum_variance)),np.arange(len(sum_variance))]),
                       'Label' : (['Explained variance']* len(sum_variance) + ['Threshold']*len(sum_variance)) })
    
    plt.figure(figsize=(10,6))
    sns.lineplot(data = df ,x = 'X',y='Explained variance',hue = 'Label',style = 'Label',marker = True)
    plt.text(len(sum_variance)-1,explained_ratio_threshold,"(%d,%0.3f)" %(len(sum_variance),sum_variance[-1]) )
    plt.grid()
    plt.xlabel("Number of principle component")
    plt.title('Explained variance ratio in different number of principle component')
    plt.savefig(output_path,dpi = 300,bbox_inches = 'tight')
    return len(sum_variance)
    
    
def centroid_base_data_augmentation(X,n_pc,n_sample_synthesized) :
    '''
    X : numpy array(n_sample,n_feature),contain origin data information
    n_pc : integer, number of principle component 
    '''
    pca = PCA(n_components=n_pc)
    pca_r = pca.fit_transform(X)
    
    n_sample = X.shape[0]
    centroid = np.sum(pca_r,axis=0) / n_sample
    distance = np.zeros(n_sample)
    '''
    p = 1 / n_sample
    p_lower = p * min(p_interval*10)/5
    p_upper = p * max(p_interval*10)/5
    p_each = (p_upper - p_lower) / (n_sample -1)
    for i in range(103) :
        sample_p[i] = p_lower + p_each * i
    '''
    #caculate probability weight for each sample
    sample_p = np.zeros(n_sample)
    #calculate distance to centroid for each data point
    for i in range(X.shape[0]) :
        distance[i] = np.linalg.norm(centroid - pca_r[i,:])
        #the closer to the centroid,the lower the probability
        sample_p[i] = i+1 / sum(np.arange(1,n_sample+1))

    synthetic_data = np.zeros([n_sample_synthesized,n_pc])
    sample_p = list(sample_p)
    
    for i in range(n_sample_synthesized) :
        choice = np.random.choice(np.argsort(distance),2,sample_p)
        r1 = pca_r[choice[0],:]
        r2 = pca_r[choice[1],:]
        synthetic_data[i,:] = (r1 + r2) / 2
    synthetic_data = pca.inverse_transform(synthetic_data)
    return synthetic_data

def data_augmentatation_pca_scatterplot(origin_m,synthetic_m,fig_output_path) :
    
    x = np.concatenate([origin_m,synthetic_m],axis=0)
    pca = PCA(n_components=3)
    syn_pca = pca.fit_transform(x)
    pca_df = pd.DataFrame({'PC1' : syn_pca[:,0],'PC2' : syn_pca[:,1],'PC3' : syn_pca[:,2],
                           'Label' : ['Origin'] *origin_m.shape[0] + ['Syntheitc'] * synthetic_m.shape[0]})

    fig,axs = plt.subplots(1,2,figsize = (16,7))
    sns.scatterplot(data = pca_df,x = 'PC1',y = 'PC2',hue = 'Label',palette='Set2',ax = axs[0])
    sns.scatterplot(data = pca_df,x = 'PC3',y = 'PC2',hue = 'Label',palette='Set2',ax=axs[1])

    axs[0].set_xlabel("PC1 (" + str(round(100*pca.explained_variance_ratio_[0],2)) + '%)')
    axs[0].set_ylabel("PC2 (" + str(round(100*pca.explained_variance_ratio_[1],2)) + '%)')
    axs[1].set_xlabel("PC3 (" + str(round(100*pca.explained_variance_ratio_[2],2)) + '%)')
    axs[1].set_ylabel("PC2 (" + str(round(100*pca.explained_variance_ratio_[1],2)) + '%)')
    plt.savefig(fig_output_path,dpi = 300,bbox_inches = 'tight')
    
    
def main() :
    parser = argparse.ArgumentParser()
    parser.add_argument("-i", "--input",help="path of expression profile")
    parser.add_argument("-r", "--reference",help="ens id for reference gene")
    parser.add_argument("-f","--fig_output",help="figure output path")
    parser.add_argument("-o","--output",help="Expression matrix with synthetic data output path")
    parser.add_argument("-n","--number",type=int,help="Number of synthetic data")
    args = parser.parse_args()
    
    exp_profile = pd.read_csv(args.input,sep='\t',index_col=0)
    hallmark = pd.read_csv(args.reference,sep='\t',header=None)
    hallmark.columns = ["EnsID"]  # type: ignore
    hallmark_gene = hallmark['EnsID'].values
    hallmark_m = exp_profile.loc[hallmark_gene,:]
    
    n_pc = plot_pc_explain_ratio(hallmark_m,args.fig_output + "principle_component_explained_ratio.png")
    synthetic_X = centroid_base_data_augmentation(hallmark_m,n_pc,args.number)
    data_augmentatation_pca_scatterplot(hallmark_m,synthetic_X,args.fig_output + "data_augmentation_pcaplot.png")
    synthetic_df = pd.DataFrame(synthetic_X,index=["Synthetic_" + str(x) for x in range(args.number)],columns=exp_profile.columns)
    pd.concat([exp_profile,synthetic_df],axis=0).to_csv(args.output,sep='\t')
    
if __name__ == '__main__' :
    main()