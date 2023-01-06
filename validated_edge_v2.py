import numpy as np
import pandas as pd
import networkx as nx
import multiprocessing as mp
from scipy import sparse
from itertools import repeat
import time
import json
import re
import argparse
import pickle

def gene_for_go(matrix,gene):
    '''
    matrix : numpy object array, columns = ['Gene stable ID', 'GO term accession', 'GO domain']
    gene : str, gene ENS ID
    '''
    if isinstance(gene,str) :
        idx = np.where(matrix[:,0] == gene,True,False)
        go = matrix[idx,1]
        return go
    else :
        return None

def create_gene2go_matrix(biomart,exp_profile,graph_dict,threads) :
    '''
    biomart : dataframe; biomart annotation of coding genes
    exp profile : dataframe; coding gene expression profile
    graph_dict : dict; record Gene ontology tree from go.obo file
    threads : int, number of threads for parallel
    '''
    df = biomart.loc[:,["Gene stable ID","GO term accession",'GO domain']]
    rm_idx = df.isnull()['GO term accession']
    biomart = df.loc[~rm_idx,:]
    gene_list = list(exp_profile.index)
    ## keep go domain and index in sparse matrix
    go2domain = dict(zip(biomart['GO term accession'],biomart['GO domain']))
    goindex = dict(zip(list(go2domain.keys()),np.arange(len(go2domain)) ))
    gene2go_matrix = sparse.dok_matrix((len(gene_list), len(goindex)),dtype=np.int) # type igonore
    ##record gene involved GO term in gene2go_dict
    m = biomart.to_numpy()
    gene_list = list(exp_profile.index)
    pool = mp.Pool(threads)
    result = pool.starmap(gene_for_go,zip(repeat(m),gene_list))
    pool.close()
    pool.join()
    gene2go_dict = dict(zip(gene_list,result))
    ##create gene to go sprase matrix
    for gene_idx,gene in enumerate(gene_list) :
        go_list = gene2go_dict[gene]
        for go in go_list :
            try :
                domain = go2domain[go]
                G = graph_dict[domain]
            except :
                continue
            try :
                anc = nx.descendants(G,go)
            except :
                continue
            try :
                index = [goindex[key] for key in anc if key in goindex.keys()] + [goindex[go]]
            except KeyError :
                index = [goindex[key] for key in anc if key in goindex.keys()]
            gene2go_matrix[gene_idx,index] = 1

    return gene2go_matrix,goindex

def get_go_level(go,go2namespace,graph_dict) :
    if go not in go2namespace.keys() :
        return -1
    field = go2namespace[go]
    root = graph_dict[field]
    return root.nodes[go]['Level']

def remove_go_by_level(gene2go_matrix,index2go,go2namespace,graph_dict,go_level_threshold=5) :
    
    go_level = np.zeros(gene2go_matrix.shape[1])
    for idx,go in enumerate(index2go) :
        go_level[idx] = get_go_level(go,go2namespace,graph_dict)
        
    idx = np.where(go_level >= go_level_threshold,1,0)
    D = sparse.diags(idx, dtype=gene2go_matrix.dtype)
    gene2go_matrix= gene2go_matrix * D
    gene2go_matrix.eliminate_zeros()
    index2go[idx == 0] = 'NA'
    
    return index2go,gene2go_matrix

from scipy.stats import fisher_exact
from math import comb

def gene_level_fisher_exact_test(gene2go_matrix,vote_result,go_index) :
    '''
    gene2go_matrix : scipy sparse matrx (csc of dok), row is gene , column is GO
    vote_result : scipy csc matrix
    '''
    matrix = gene2go_matrix.getcol(go_index)
    if matrix.count_nonzero() == 0 :
        p = 999
        return p
    intersection,_ = matrix.nonzero()
    vote,_ = vote_result.nonzero()
    g1 = len(np.intersect1d(intersection,vote))
    g2 = len(intersection) - g1
    g3 = len(vote) - g1
    g4 = gene2go_matrix.shape[0] - len(intersection - g3)

    fisher_table = np.array([[g1,g2],[g3,g4]])
    try :
        oddsr, p = fisher_exact(fisher_table, alternative='two-sided')
    except ValueError :
        print("Error occur in fisher exact test stage !")
        p = 999

    return p

def gene_based_fisher_exact_test(gene2go_matrix,vote_matrix) :
    '''
    gene2go_matrix : scipy csr matrix; contain go & gene interaction 
    vote : np.array ; Ensemble learning vote result of each gene
    '''
    vote_csc = sparse.csc_matrix(vote_matrix)
    vote_csc.eliminate_zeros()
    
    n = gene2go_matrix.shape[1]
    gene_fisher_pv = np.zeros(n)
    s = time.time()
    for i in range(n) :
        gene_fisher_pv[i] = gene_level_fisher_exact_test(gene2go_matrix,vote_csc,i)
        if i % 1000 == 0 :
            print('Number of %d GO fisher exact test is processed !' % i)
    e = time.time()
    delta = e -s
    print("Execution time of gene based fisher exact test is %0.2f seconds" % delta)

    pass_idx = np.where(gene_fisher_pv < 0.05,1,0)
    D = sparse.diags(pass_idx, dtype=gene2go_matrix.dtype)
    gene2go_matrix = gene2go_matrix * D
    gene2go_matrix.eliminate_zeros()
    
    return gene2go_matrix,pass_idx

def corr_2_sparse_matrix(corr_matrix,threshold=0.4) :
    corr_matrix = np.tril(corr_matrix)
    idx = corr_matrix >= threshold  # type: ignore
    corr_matrix[~idx] = 0  # type: ignore
    np.fill_diagonal(corr_matrix,0)
    corr_csc = sparse.csc_matrix(corr_matrix)
    corr_csc.eliminate_zeros()
    return corr_csc

def sum_csc_matrix(matrix,index) :
    tmp = matrix[:,index]
    tmp = tmp.tocsr()
    nonzero = 0
    for i in index :
        nonzero +=  tmp.getrow(i).count_nonzero()
    del tmp
    return nonzero

def edge_fisher_exact_test(gene2go_matrix,corr_matrix,go_index) :
    
    matrix = gene2go_matrix.getcol(go_index)
    if matrix.count_nonzero() == 0 :
        p = 999 
        return p
    else :
        #GO involved gene set
        gene_idx,_ = matrix.nonzero()
        num_gene_involved = len(gene_idx)
        m = corr_matrix[gene_idx][:,gene_idx]
        e1 = m.count_nonzero()
        #e1 = sum_csc_matrix(corr_matrix,gene_idx)
        #e1 + e2
        total_intersection = comb(num_gene_involved,2)
        e2 = total_intersection - e1
        ### GO with intersection
        #e1 + e3
        total_corr = corr_matrix.count_nonzero()
        e3 = total_corr - e1
        n_gene = corr_matrix.shape[0]
        e4 = comb(n_gene,2) - comb(num_gene_involved,2) - e3

        fisher_table = np.array([[e1,e2],[e3,e4]])
        try :
            _, p = fisher_exact(fisher_table, alternative='two-sided')
        except ValueError :
            print("Error occur in fisher exact test stage !")
            p = 999

        return p
    
def create_go2gene_dict(gene2go_matrix,go_list,biomart) :
    '''
    gene2go_matrix : scipy csr matrix; row is gene, column is GO term
    go_list : numpy array, go term of gene2go_matrix
    biomart : dataframe; biomart annotation information 
    '''
    go2gene_dict = dict()
    idx = np.where(np.sum(gene2go_matrix,axis=0)[0,:] > 0)[1]  # type: ignore
    validated_go = go_list[idx]
    ###
    missing = 0
    n = 0
    for g in validated_go :
        n += 1
        idx = np.where(biomart['GO term accession'].values == g,True,False)
        if sum(idx) > 0 :
            genes = biomart['Gene stable ID'][idx]
            genes = np.unique(genes)
            go2gene_dict[g] = genes
        else :
            missing += 1
        if n % 200 == 0 :
            print(n)
    return go2gene_dict


def main() :
    
    parser = argparse.ArgumentParser()
    parser.add_argument("-i", "--input",help="path of lihc expression profile")
    parser.add_argument("-b", "--biomart",help="path of biomart information")
    parser.add_argument("-g", "--graph",help="path of go graph object")
    parser.add_argument("-v","--vote",help="path of vote result")
    parser.add_argument("-c", "--correlation_threshold",help="threshold of correlation between 2 nodes")
    parser.add_argument("-t", "--threads",type=int,help="number of threads to use")
    parser.add_argument("-o","--output_path",help = 'path of validated go information')
    parser.add_argument("-n","--name",help='name of project')
    args = parser.parse_args()
    
    exp_profile = pd.read_csv(args.input,sep='\t',index_col = 0)
    gene_idx = [bool(re.search("Synthetic",x)) == False for x in exp_profile.index]
    exp_profile = exp_profile.loc[gene_idx,:]
    biomart = pd.read_csv(args.biomart,sep='\t')
    #import go graph
    graph_dict = dict()
    #graph_path = '/home/bruce1996/data/GO/networkx/'
    go_field = ['biological_process','cellular_component','molecular_function']
    for field in go_field :
        file = args.graph + field + '.go'
        with open(file , 'rb') as f:  # notice the r instead of w
            graph_dict[field] = pickle.load(f)
            f.close()
    with open("/home/bruce1996/data/GO/go2namespace.json",'rb') as f :
        go2namespace = json.load(f)
        f.close()
    # remove go information by biomart evidence code
    experiment_evidence_code = ['EXP','IDA','IPI','IMP','IGI','IEP','HTP','HDA','HMP','HGI','HEP']
    idx = [x in experiment_evidence_code for x in biomart['GO term evidence code']]
    biomart = biomart.loc[idx,:]
            
    print("Gene2go stage processing!")
    gene2go_matrix,goindex = create_gene2go_matrix(biomart,exp_profile,graph_dict,threads = args.threads)
    index2go = np.array(list(goindex.keys()))
    gene2go_matrix = gene2go_matrix.tocsc()
    print("Calculate correlation matrix !")
    corr_m = exp_profile.T.corr(method='spearman').to_numpy()
    ## remove go term which level less than threshold
    index2go,gene2go_matrix = remove_go_by_level(gene2go_matrix,index2go,go2namespace,graph_dict)
    ## import vote file
    vote = pd.read_csv(args.vote,sep='\t',index_col=0)
    vote = vote.loc[exp_profile.index,:]
    vote_csc = sparse.csc_matrix(vote.to_numpy())
    vote_csc.eliminate_zeros()
    ## gene level fisher exact test
    n = gene2go_matrix.shape[1]
    gene_fisher_pv = np.zeros(n)
    s = time.time()
    for i in range(n) :
        gene_fisher_pv[i] = gene_level_fisher_exact_test(gene2go_matrix,vote_csc,i)
        if i % 1000 == 0 :
            print('Number of %d GO fisher exact test is processed !' % i)
    e = time.time()
    delta = e -s
    print("Execution time of gene based fisher exact test is %0.2f seconds" % delta)

    idx = np.where(gene_fisher_pv < 0.05,1,0)
    D = sparse.diags(idx, dtype=gene2go_matrix.dtype)
    gene2go_matrix_pass_gene = gene2go_matrix * D
    gene2go_matrix_pass_gene.eliminate_zeros()  # type: ignore
    ## edge level fisher exact test
    corr_csc = corr_2_sparse_matrix(corr_m)
    n = gene2go_matrix_pass_gene.shape[1]  # type: ignore
    edge_fisher_pv = np.zeros(n)
    s = time.time()
    for i in range(n) :
        edge_fisher_pv[i] = edge_fisher_exact_test(gene2go_matrix_pass_gene,corr_csc,i)
        if i % 500 == 0 :
            print('Number of %d GO fisher exact test is processed !' % i)
    e = time.time()
    delta = e -s
    print("Execution time of gene based fisher exact test is %0.2f second" % delta)

    idx = np.where(edge_fisher_pv < 0.05,1,0)
    D = sparse.diags(idx, dtype=gene2go_matrix_pass_gene.dtype)  # type: ignore
    gene2go_matrix_validated = gene2go_matrix_pass_gene * D
    gene2go_matrix_validated.eliminate_zeros()

    sparse.save_npz(args.output_path + args.name + '_validated_gene2go_matrix.npz',gene2go_matrix_validated)
    
    gene_in_go = create_go2gene_dict(gene2go_matrix_validated,index2go,biomart)
    with open(args.output_path + args.name +"_validated_go_2_gene_dict.pkl",'wb') as f :
        pickle.dump(gene_in_go,f)
    f.close()
    
if __name__ == '__main__' :
    main()
        