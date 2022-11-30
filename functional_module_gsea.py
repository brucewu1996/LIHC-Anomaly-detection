import numpy as np
import pandas as pd
import pickle,time
import argparse
from scipy import sparse
from utils.gsea import enrich_score,permutation,gene_set_enrichment_analysis


class gene_set_enrichment(gene_set_enrichment_analysis) :
    
    def gsea_with_weight(self,weight_array) :
        '''
        weight_array : numpy array, contain only 1/0
        '''
        component_list = np.array(list(self.cluster_component.keys()))
        n_component = len(self.cluster_component.keys())
        include_idx = np.where(weight_array != 0,True,False)
        self.ranking = np.array(self.ranking)[include_idx]

        for idx,cluster in enumerate(component_list) :
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
    '''
    vote_csc = sparse.csc_matrix(vote.to_numpy())
    vote_csc.eliminate_zeros()
    '''
    ### vote level GSEA 
    print("GSEA analysis by vote of functional module is processing !!")
    vote_start = time.time()
    vote_ranking = list(vote.sort_values(by='Vote',ascending=False).index)
    vote_weight = np.where(vote.sort_values(by='Vote',ascending=False)['Vote'].values > 0 ,1,0)  # type: ignore
    vote_tsea = gene_set_enrichment(go_genedict,vote_ranking)
    vote_tsea.gsea_with_weight(vote_weight)
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