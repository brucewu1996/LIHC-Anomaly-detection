import numpy as np
from numba import jit

@jit(nopython=True)
def enrich_score(taxa_array) : 
    '''
    pathway : name of pathway, ex : DTDPRHAMSYN-PWY: dTDP-&beta;-L-rhamnose biosynthesis
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

@jit(nopython=True)
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

class taxonomy_set_enrichment_analysis :
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
        
    def tsea(self) :
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
        
