import numpy as np
import pandas as pd
import pickle,time
import argparse
from scipy import sparse
from utils.gsea import permutation,gene_set_enrichment_analysis,matrix_2_corr_sparse_matrix,remove_uncorr_component


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
    parser.add_argument("-g","--go",help="GO term 2 geneset dict after fisher exact test")
    parser.add_argument("-v","--vote",help="path of vote result")
    parser.add_argument("-d","--deseq",help="prefix of vote result")
    parser.add_argument("-e","--exp_profile",help='path of expression profile')
    parser.add_argument("-o","--output_path",type=str,help = 'path of ensemble model output')
    parser.add_argument("-p","--prefix",help = 'prefix of figure')
    parser.add_argument("-w","--weight",type=str,help='GSEA with vote only gene or Not')
    args = parser.parse_args()  
    
    ### load GSEA go to gene dictionary
    with open(args.go,'rb') as f :
        go_genedict = pickle.load(f)
    exp_profile = pd.read_csv(args.exp_profile,sep='\t',index_col=0)
    #remove gene without any correlation > 0.4 with other genes
    remove_uncorr_component(go_genedict,exp_profile)
    ###load GSEA sorting file (vote & fold change)
    ### vote level GSEA 
    print("GSEA analysis by vote of functional module is processing !!")
    vote = pd.read_csv(args.vote,sep='\t',index_col=0)
    vote_start = time.time()
    if args.weight == "True" :
        vote_ranking = list(vote.sort_values(by='Vote',ascending=False).index)
        vote_weight = np.where(vote.sort_values(by='Vote',ascending=False)['Vote'].values > 0 ,1,0)  # type: ignore
        print("Number of voted gene : %d" % sum(vote_weight))
        vote_tsea = gene_set_enrichment(go_genedict,vote_ranking)
        vote_tsea.gsea_with_weight(vote_weight)
    else :
        vote_ranking = list(vote.sort_values(by='Vote',ascending=False).index)
        vote_tsea = gene_set_enrichment(go_genedict,vote_ranking)
        vote_tsea.gsea()
    print(sum(vote_tsea.pesudo_f))
    df = vote_tsea.validate_component(pesudo_f_threshold=0.3)
    pass_vote_go = list(df.index)
    print("Number of GO pass vote-driven GSEA analysis : %d" % len(pass_vote_go))
    pass_vote_go_gene_dict = dict()
    for go in pass_vote_go :
        pass_vote_go_gene_dict[go] = go_genedict[go]
    vote_end = time.time()
    print("Execution time of GSEA analysis by vote number is : %0.2f seconds" % (vote_end - vote_start))
    #activate functional module
    #fold change
    print("GSEA analysis by fold-change of activate functional module is processing !!")
    deseq_df = pd.read_csv(args.deseq,sep='\t',index_col=0)
    fc_activated_start = time.time()
    fc_ranking = list(deseq_df.sort_values(by='log2FoldChange',ascending=False).index)
    foldchange_gsea_activated = gene_set_enrichment_analysis(pass_vote_go_gene_dict,fc_ranking)
    foldchange_gsea_activated.gsea()
    validated_df = foldchange_gsea_activated.validate_component()
    print("Number of GO pass fold-change-driven GSEA analysis (activated) : %d" % validated_df.shape[0])
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
    print("Number of GO pass fold-change-driven GSEA analysis (activated) : %d" % validated_df.shape[0])
    validated_df.to_csv(args.output_path + args.prefix + "_inactivated_functional_module.txt",sep='\t')
    fc_inactivated_end = time.time()
    print("Execution time of GSEA analysis of activated functional module is : %0.2f seconds" % (fc_inactivated_end - fc_inactivated_start))
    print("Execution time of whole GSEA analysis of functional module is : %0.2f seconds" % (fc_inactivated_end - vote_start))
    
if __name__ == '__main__' :
    main()