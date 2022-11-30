import numpy as np
import pandas as pd
import argparse,pickle,os,time
from scipy import sparse
from numba import jit


@jit(nopython=True)  # type: ignore
def enrich_score(taxa_array) : 
    '''
    taxonomy_list : taxonomy list with "order"
    '''
    #calculate enrichment score unit, enrichment score upper bond = 1, lower bond = -1
    es = 0
    max_es = 0
    if sum(taxa_array) == 0 : 
        return max_es
    elif (len(taxa_array) - sum(taxa_array)) == 0 :
        max_es = 1
        return max_es
    else :
        es_pos = 1/sum(taxa_array)
        es_n = 1/(len(taxa_array) - sum(taxa_array))
        #calculate enrichment score by "ordered" taxonomy list
        for t in taxa_array :
            if t  == 1 :
                es += es_pos
            elif t == 0 :
                es -= es_n
            else :
                print('Value in taxa array must be 0 or 1')
                break
            if es > max_es :
                max_es = es
    return max_es

@jit(nopython=True)  # type: ignore
def permutation(taxa_array,n_times = 1000) :
    
    origin_es = enrich_score(taxa_array)
    permutation_es = np.zeros(n_times)
    for i in range(n_times) :
        np.random.shuffle(taxa_array)
        permutation_es[i] = enrich_score(taxa_array)
    es_above_origin = sum(permutation_es > origin_es)
    if es_above_origin == 0 :
        pesudo_f = 1/ (n_times * 10)
    else :
        pesudo_f = es_above_origin / n_times
    return pesudo_f,origin_es

class gene_set_enrichment_analysis :
    def __init__(self,cluster_component,ranking,permutaion= 1000) :
        '''
        cluster_component : dict; key equal to cluster name,value equal to cluster member
        ranking : list; sorted component 
        permutation : int; number of permutation 
        '''
        self.cluster_component = cluster_component
        self.ranking = ranking
        self.permutation = permutaion
        self.es_score = np.zeros(len(cluster_component))
        self.pesudo_f = np.zeros(len(cluster_component))
        
    def gsea(self) :
        n_component = len(self.cluster_component.keys())
        for idx,cluster in enumerate(self.cluster_component.keys()) :
            target = self.cluster_component[cluster]
            es_array = np.zeros(len(self.ranking))
            for r_idx,r in enumerate(self.ranking) :
                if r in target :
                    es_array[r_idx] = 1
            if sum(es_array) > 0 :
                self.pesudo_f[idx],self.es_score[idx] = permutation(es_array,n_times=self.permutation)
            else :
                self.pesudo_f[idx],self.es_score[idx] = (1,0)
            
            if idx % 50 == 0 and idx > 0:
                print("Number of components are processed : %d / %d" % (idx,n_component))
                
    def validate_component(self,pesudo_f_threshold=0.05) :
        
        component_list = list(self.cluster_component.keys())
        df = pd.DataFrame({'Pesudo-F' : self.pesudo_f},index = component_list)
        idx = np.where(df['Pesudo-F'] < pesudo_f_threshold,True,False)
        validated_df = df.loc[idx,:].sort_values(by='Pesudo-F')
        return validated_df
                
def main() :
    parser = argparse.ArgumentParser()
    parser.add_argument("-g", "--go",help="GO term 2 geneset dict after fisher exact test")
    parser.add_argument("-v", "--vote",help="path of vote result")
    parser.add_argument("-d", "--deseq",help="prefix of vote result")
    parser.add_argument("-o","--output_path",type=str,help = 'path of ensemble model output')
    parser.add_argument("-p","--prefix",help = 'prefix of figure')
    args = parser.parse_args()  
    
    # load validated gene2go matrix
    #gene2go_matrix_validated = sparse.load_npz('/home/bruce1996/data/LIHC_anomaly_detection/functional_profiling/candidate_functional_module/validated_gene2go_matrix.npz')
    # load validate gene2go dict
    #/home/bruce1996/data/LIHC_anomaly_detection/functional_profiling/candidate_functional_module/validated_go2gene_dict.pkl'
    with open(args.go,'rb') as f :
        go_genedict = pickle.load(f)
    ###load GSEA sorting file (vote & fold change)
    # vote
    vote = pd.read_csv(args.vote,index_col=0)
    vote_csc = sparse.csc_matrix(vote.to_numpy())
    vote_csc.eliminate_zeros()
    ### vote level GSEA 
    print("GSEA analysis by vote of functional module is processing !!")
    vote_start = time.time()
    vote_ranking = list(vote.sort_values(by='Vote',ascending=False).index)
    vote_tsea = gene_set_enrichment_analysis(go_genedict,vote_ranking)
    vote_tsea.gsea()
    df = vote_tsea.validate_component()
    pass_vote_go = list(df.index)
    pass_vote_go_gene_dict = dict()
    for go in pass_vote_go :
        pass_vote_go_gene_dict[go] = go_genedict[go]
    vote_end = time.time()
    print("Execution time of GSEA analysis by vote number is : %0.2f seconds" % (vote_end - vote_start))
    # fold change
    deseq_df = pd.read_csv(args.deseq,sep='\t',index_col=0)
    #activate functional module
    print("GSEA analysis by fold-change of activate functional module is processing !!")
    fc_activated_start = time.time()
    fc_ranking = list(deseq_df.sort_values(by='log2FoldChange',ascending=False).index)
    foldchange_gsea_activated = gene_set_enrichment_analysis(pass_vote_go_gene_dict,fc_ranking)
    foldchange_gsea_activated.gsea()
    validated_df = foldchange_gsea_activated.validate_component()
    validated_df.to_csv(args.output_path + args.prefix + "_activated_functional_module.txt",sep='\t')
    fc_activated_end = time.time()
    print("Execution time of GSEA analysis of activated functional module is : %0.2f seconds" % (fc_activated_end - fc_activated_start))
    #suppress functional module
    print("GSEA analysis by fold-change of inactivate functional module is processing !!")
    fc_inactivated_start = time.time()
    fc_ranking = list(deseq_df.sort_values(by='log2FoldChange').index)
    foldchange_gsea_inactivated = gene_set_enrichment_analysis(pass_vote_go_gene_dict,fc_ranking)
    foldchange_gsea_inactivated.gsea()
    validated_df = foldchange_gsea_inactivated.validate_component()
    validated_df.to_csv(args.output_path + args.prefix + "_inactivated_functional_module.txt",sep='\t')
    fc_inactivated_end = time.time()
    print("Execution time of GSEA analysis of activated functional module is : %0.2f seconds" % (fc_inactivated_end - fc_inactivated_start))
    print("Execution time of whole GSEA analysis of functional module is : %0.2f seconds" % (fc_inactivated_end - vote_start))
    
if __name__ == '__main__' :
    main()
        
