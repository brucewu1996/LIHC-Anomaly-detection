import sys,pickle,json,time
sys.path.append("/home/bruce1996/repo/LIHC_anomaly_detection/")
from utils.gsea import *
from utils.go2gene_dict import *
from numba import jit
import numpy as np
import pandas as pd
import networkx as nx
import argparse

def format_gsea_result(gsea_result,go_graph_dict,go2namespace) :

    component_namespace,component_level,component_description = [],[],[]
    for go in gsea_result.index :
        namespace = go2namespace[go]
        go_level = go_graph_dict[namespace].nodes[go]['Level']
        des = go_graph_dict[namespace].nodes[go]['name']
        component_namespace.append(namespace)
        component_description.append(des)
        component_level.append(go_level)
    gsea_result['GO_level'] = component_level
    gsea_result['Description'] = component_description
    gsea_result['namespace'] = component_namespace

    column_order = ['GO_level','#Nodes','namespace','Description','Pesudo-F','Origin_ES','ES_median','ES_std']
    gsea_result = gsea_result.loc[:,column_order]
    return gsea_result

def main() :

    parser = argparse.ArgumentParser()
    parser.add_argument("-d", "--deseq",help="path of deseq result")
    parser.add_argument("-p","--output_prefix",help="output file prefix")
    parser.add_argument("-o","--output",help="Expression matrix with synthetic data output path")
    args = parser.parse_args()
    ### load go graph
    namespace = ['biological_process','cellular_component','molecular_function']
    graph_d = dict()
    for domain in namespace :
        with open("/home/bruce1996/data/GO/networkx/%s.go" % domain,'rb') as f :
            g = pickle.load(f)
        f.close()
        graph_d[domain] = g
    ### load go2gene dict
    with open("/home/bruce1996/data/GO/networkx/go2gene_dict.pkl",'rb') as f :
        go2gene_dict = pickle.load(f)
    f.close()
    with open("/home/bruce1996/data/GO/go2namespace.json",'rb') as f :
        go2namespace = json.load(f)
    f.close()
    ###
    deseq_res = pd.read_csv(args.deseq,sep='\t',index_col=0)
    ascending_ranking = list(deseq_res.sort_values(by='log2FoldChange').index)
    descending_ranking = ascending_ranking[::-1]
    ###activate gsea
    start_gsea = time.time()
    gsea = gene_set_enrichment_analysis(go2gene_dict,descending_ranking)
    gsea.gsea()
    res = gsea.validate_component()
    gsea_res = format_gsea_result(res,graph_d,go2namespace)
    gsea_res.to_csv(args.output + "%s_activated_gsea_summary.txt" % args.output_prefix,sep='\t')
    end_gsea = time.time()
    exe_time = end_gsea - start_gsea
    print("Execution time of GSEA is %d : %0.2f " % (exe_time // 60, exe_time % 60))
    ###activate gsea
    start_gsea = time.time()
    gsea = gene_set_enrichment_analysis(go2gene_dict,ascending_ranking)
    gsea.gsea()
    res = gsea.validate_component()
    gsea_res = format_gsea_result(res,graph_d,go2namespace)
    gsea_res.to_csv(args.output + "%s_suppressed_gsea_summary.txt" % args.output_prefix,sep='\t')
    end_gsea = time.time()
    exe_time = end_gsea - start_gsea
    print("Execution time of GSEA is %d mins %0.2f seconds " % (exe_time // 60, exe_time % 60))

if __name__ == '__main__' :
    main()