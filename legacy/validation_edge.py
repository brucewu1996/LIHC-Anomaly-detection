import numpy as np
import pandas as pd
import networkx as nx
import multiprocessing as mp
from scipy import sparse
from itertools import repeat,combinations
from functools import partial
import tracemalloc
import time
import os
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

def candidate_edge(gene_go_matrix,edge) :
    '''
    cadidate edge is gene pair in same GO term which level under 5
    '''
    flag = False
    gene1,gene2 = edge
    overlap = gene_go_matrix[gene1,:] + gene_go_matrix[gene2,:]
    #find non-zero element in sparse matrix in gene 1/2 row
    x,y,val = sparse.find(overlap)
    # val == 2 means both gene in this GO term is non-zero
    idx = np.where(val == 2)[0]
    if len(idx) > 0 :
        flag =True
    return flag,idx

def split_edge_list(edges,n_block) :
    '''
    split edge list into many pieces for parallel computing
    '''
    sample_per_part = round(len(edges) / n_block)
    edge_list = list()
    for i in range(n_block):
        start = i * sample_per_part
        end = (i+1) * sample_per_part
        if end > len(edges) :
            end = len(edges) - 1
        partial_edge = edges[start:end]
        edge_list.append(partial_edge)
    return edge_list

def edge_list_with_corr(edge_list,corr_matrix) :
    edge_with_corr = list()
    for i in range(len(edge_list)) :
        edge_dict = dict()
        for e in edge_list[i] :
            edge_dict[e] = corr_matrix[e[0],e[1]]
        edge_with_corr.append(edge_dict)
    return edge_with_corr

def validate_edge(edges,gene2go_matrix,index2go,corr_threshold = 0.4) :
    '''
    edges : list, produce by itertools combination, key is edges (index1,index2),val is correlation (ex : 0.4xxxx)
    corr_matrix : numpy array, produce by df.corr()
    gene2go_matrix : scipy sparse matrix, record gene involved GO term, row is gene col is go
    index2go : numpy array, record repensented GO term in gene2go_matrix index (column)
    '''
    intersection_go = dict()
    print(f'validation edge : {mp.current_process()=}')
    #dot_matrix = dok_matrix((corr_matrix.shape[0],corr_matrix.shape[1]),dtype=np.int8)
    for edge,value in edges.items() :
        x = edge[0]
        y = edge[1]
        #filter by correlation
        if value > corr_threshold :
            #filter by intersection go term
            flag,go_pos = candidate_edge(gene2go_matrix,edge)
            if flag :
                intersection_go[edge] = index2go[go_pos]
        else :
            continue
    return intersection_go

def save_edge_ancestor(result,output_path) :
    for i in range(len(result)) :
        file = output_path + 'validation_edge_' + str(i) + '.pkl'
        with open(file, 'wb') as fp:
            pickle.dump(result[i], fp)

def save_edge2graph(result,gene_list,corr_matrix,output_path) :
    G = nx.DiGraph()
    for edge_dict in result :
        graph_dict = dict()
        for edge in edge_dict.keys() :
            n1 = gene_list[edge[0]]
            n2 = gene_list[edge[1]]
            corr = corr_matrix[edge[0],edge[1]]
            graph_dict[(n1,n2)] = corr
        G.add_edges_from(graph_dict)
        nx.set_edge_attributes(G,graph_dict,'Correlation')
    file = output_path + 'validation_gene_correlation_network.pkl'
    with open(file,'wb') as fp :
        pickle.dump(G, fp)

def main() :

    parser = argparse.ArgumentParser()
    parser.add_argument("-i", "--input",help="path of lihc expression profile")
    parser.add_argument("-b", "--biomart",help="path of biomart information")
    parser.add_argument("-g", "--graph",help="path of go graph object")
    parser.add_argument("-c", "--correlation_threshold",help="threshold of correlation between 2 nodes")
    parser.add_argument("-t", "--threads",type=int,help="number of threads to use")
    parser.add_argument("-o","--output_path",help = 'path of ensemble model output')
    parser.add_argument("-s","--stand",type=bool,default=False,help = "Standardziation of not,if standardzied choose False")
    args = parser.parse_args()
    #load essenial files
    exp_profile = pd.read_csv(args.input,sep='\t',index_col = 0)
    biomart = pd.read_csv(args.biomart,sep='\t')
    #import go graph
    graph_dict = dict()
    go_field = ['biological_process','cellular_component','molecular_function']
    for field in go_field :
        file = args.graph + field + '.go'
        with open(file , 'rb') as f:  # notice the r instead of w
            graph_dict[field] = pickle.load(f)
    #create a gene2go matrix
    print("Gene2go stage processing!")
    gene_list = list(exp_profile.index)
    gene2go_matrix,goindex = create_gene2go_matrix(biomart,exp_profile,graph_dict,threads = args.threads)
    index2go = np.array(list(goindex.keys()))
    print("Calculate correlation matrix !")
    index = np.arange(exp_profile.shape[0],dtype=int)
    edges = list(combinations(index,2))
    corr_m = exp_profile.T.corr(method='spearman').to_numpy()
    #split edges into 1000 partition with correlation value
    print("Split edge list!")
    edge_list = split_edge_list(edges,1000)
    edge_with_corr_list = edge_list_with_corr(edge_list,corr_m)
    gene2go_csr = gene2go_matrix.tocsr()
    tracemalloc.start()
    start_time = time.time()
    #multiprocessing
    pool = mp.Pool(args.threads)
    corr_threshold = float(args.correlation_threshold)
    result = pool.map(partial(validate_edge,gene2go_matrix=gene2go_csr,index2go=index2go,corr_threshold=corr_threshold),edge_with_corr_list)
    pool.close()
    pool.join()
    #memory record
    current, peak = tracemalloc.get_traced_memory()
    print(f"Current memory usage {current/1e6}MB; Peak: {peak/1e6}MB")
    print(f'Time elapsed: {time.time()-start_time:.2f}s')
    tracemalloc.stop()
    #make sure all children process be terminated
    active = mp.active_children()
    for c in active :
        c.terminate()

    if os.path.exists(args.output_path) == False :
        os.mkdir(args.output_path)
    save_edge_ancestor(result,args.output_path)
    save_edge2graph(result,gene_list,corr_m,args.output_path)

if __name__ == '__main__' :
    main()