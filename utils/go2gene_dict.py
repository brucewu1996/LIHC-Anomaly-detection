import numpy as np
import pandas as pd
import networkx as nx
import multiprocessing as mp
from scipy import sparse
from itertools import repeat
from numba import jit

class gene_go_matrix :

    def __init__(self,rowname,colname,matrix) :
        self.matrix = matrix
        self.index = rowname
        self.columns = colname

def create_hgsc2ens_converter(ensembel_info,ens_col = 'gene_id',hgsc_col = 'gene_name') :
    '''
    ensembel_info : str; path of ensembel information
    '''
    info = pd.read_csv(ensembel_info,sep='\t',index_col=0)
    hgsc2ens = dict(zip(info[hgsc_col],info[ens_col]))
    return hgsc2ens

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

def create_gene2go_matrix(biomart,graph_dict,threads) :
    '''
    biomart : dataframe; biomart annotation of coding genes
    graph_dict : dict; record Gene ontology tree from go.obo file
    threads : int, number of threads for parallel
    '''
    gene_list = list(set(biomart["Gene stable ID"]))
    # remove go information by biomart evidence code
    experiment_evidence_code = ['EXP','IDA','IPI','IMP','IGI','IEP','HTP','HDA','HMP','HGI','HEP']
    idx = [x in experiment_evidence_code for x in biomart['GO term evidence code']]
    biomart = biomart.loc[idx,:]
    idx = np.where((biomart['GO domain'] != np.NaN) & (biomart['GO term accession'].isnull() == False),True,False)
    biomart = biomart.loc[idx,["Gene stable ID","GO term accession",'GO domain']]   
    ## keep go domain and index in sparse matrix
    go2domain = dict(zip(biomart['GO term accession'],biomart['GO domain']))
    goindex = dict(zip(list(go2domain.keys()),np.arange(len(go2domain)) ))
    gene2go_matrix = sparse.dok_matrix((len(gene_list), len(goindex)),dtype=np.int)  # type: ignore
    ##record gene involved GO term in gene2go_dict
    m = biomart.to_numpy()
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
                anc = nx.descendants(G,go)  # type: ignore
            except :
                continue
            try :
                index = [goindex[key] for key in anc if key in goindex.keys()] + [goindex[go]]
            except KeyError :
                index = [goindex[key] for key in anc if key in goindex.keys()]
            gene2go_matrix[gene_idx,index] = 1
    # assign gene2go_matrix into class 
    genego_matrix = gene_go_matrix(rowname=np.array(gene_list),colname=np.array(list(goindex.keys())),matrix = gene2go_matrix)

    return genego_matrix

def go2namespace(biomart) :
    
    experiment_evidence_code = ['EXP','IDA','IPI','IMP','IGI','IEP','HTP','HDA','HMP','HGI','HEP']
    idx = [x in experiment_evidence_code for x in biomart['GO term evidence code']]
    biomart = biomart.loc[idx,:]
    idx = np.where((biomart['GO domain'] != np.NaN) & (biomart['GO term accession'].isnull() == False),True,False)
    biomart = biomart.loc[idx,["Gene stable ID","GO term accession",'GO domain']]
    ## keep go domain and index in sparse matrix
    go2domain = dict(zip(biomart['GO term accession'],biomart['GO domain']))
    return go2domain

def get_go_level(go,go2namespace,graph_dict) :
    if go not in go2namespace.keys() :
        return -1
    field = go2namespace[go]
    root = graph_dict[field]
    return root.nodes[go]['Level']

def remove_go_by_level(gene_go_matrix,go2namespace,graph_dict,go_level_threshold=5) :
    '''
    gene2go_matrix : class gene_go_matrix, gene_go_matrix.matrix is spipy.csr matrix, 
    gene_go_matrix.index is gene name of matrix, gene_go_matrix.columns is go term of matrix
    go2namespace : dict, key is GO term for root of each go domain 
    graph_dict : dict, record tree of each GO domain
    go_level_threshold, int, go under this integer will excluded
    '''
    
    gene2go_matrix = gene_go_matrix.matrix
    go_level = np.zeros(gene2go_matrix.shape[1])
    for idx,go in enumerate(gene_go_matrix.columns) :
        go_level[idx] = get_go_level(go,go2namespace,graph_dict)
        
    idx = np.where(go_level >= go_level_threshold,1,0)
    D = sparse.diags(idx, dtype=gene2go_matrix.dtype)
    gene2go_matrix= gene2go_matrix * D
    gene2go_matrix.eliminate_zeros()
    gene_go_matrix.columns[idx == 0] = 'NA'
    gene_go_matrix.matrix = gene2go_matrix
    
    return gene_go_matrix

def go2dict(go,biomart) :
    go2gene_dict = {}
    idx = np.where(biomart['GO term accession'].values == go,True,False)
    if sum(idx) > 0 :
        genes = biomart['Gene stable ID'][idx]
        genes = np.unique(genes)
        go2gene_dict[go] = list(genes)
    return go2gene_dict

def create_go2gene_dict(gene_go_matrix,biomart,threads=24) :
    '''
    gene2go_matrix : go_gene_matrix class, go_gene_matrix.matrix is scipy csr matrix; 
    row is gene, column is GO term
    biomart : dataframe; biomart annotation information 
    '''
    go2gene_dict = {}
    idx = np.where(np.sum(gene_go_matrix.matrix,axis=0)[0,:] > 0)[1]  # type: ignore
    validated_go = gene_go_matrix.columns[idx]
    ###
    pool = mp.Pool(threads)
    result = pool.starmap(go2dict,zip(validated_go,repeat(biomart)))
    pool.close()
    pool.join()
    #merge multiple processing result
    for d in result:
        for k, v in d.items():  # d.items() in Python 3+
            go2gene_dict[k] = v

    return go2gene_dict
