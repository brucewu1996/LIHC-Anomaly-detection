import numpy as np
import pandas as pd
import os,pickle
import networkx as nx

def obo2dict(file,verbose=True) :
    '''
    save single go term obo file information as dictionary
    file : str; path of obo file
    '''
    #remove attr def,comment and synonym
    single_value = ['id','name','namespace','alt_id']
    multi_value = ['subset','xref','is_a','relationship']
    attr_dict = dict()
    attr_dict['edge'] = dict()

    f = open(file, 'r')
    text = f.read()
    #check go term is valided
    if 'is_obsolete' in text :
        if verbose :
            print("This go term is depercated !")
        return dict()

    for line in text.splitlines():
        if not ': ' in line:
            continue
        else :
            key, val = line.split(': ',1)
            if key in single_value :
                attr_dict[key] = val
            elif key in multi_value :
                #if key == is_a | relationship,change val format
                if key != 'is_a' and key != 'relationship' :
                    if key not in attr_dict.keys():
                        attr_dict[key] = [val]
                    else:
                        attr_dict[key].append(val)
                else :
                    if 'id' not in attr_dict.keys():
                        continue
                    elif key == 'is_a' :
                        val = val.split(' ! ',1)[0]
                        id = attr_dict['id']
                        edge_key = (id,val)
                        attr_dict['edge'][edge_key] = 'is_a'
                    elif key == "relationship" and "part_of" in val :
                        val = val.split(' ! ',1)[0].split()[1]
                        id = attr_dict['id']
                        edge_key = (id,val)
                        attr_dict['edge'][edge_key] = 'part_of'
                    elif key == "relationship" and "part_of" not in val :
                        continue

    return attr_dict

def obo2graph(obo_path) :

    name_dict = dict()
    namespace_dict = dict()
    edge = dict()
    edge_dict = {'name' : name_dict , 'namespace' : namespace_dict,'edge' : edge}
    # variable for graph
    bp = nx.DiGraph()
    cc = nx.DiGraph()
    mf = nx.DiGraph()
    graph_dict = {'biological_process' : bp,'cellular_component' : cc ,'molecular_function' : mf}

    #convert obo to dict
    go_list = os.listdir(obo_path)
    for go in go_list :
        f = obo_path + go
        obo_dict = obo2dict(f,verbose=False)
        #keep GO Term which belong to bp,cc and mf
        try :
            if obo_dict['namespace'] not in graph_dict.keys() :
                continue
        except KeyError :
                continue
        #save attr to edgedict
        if 'id' in obo_dict.keys() :
            id = obo_dict['id']
            for key in obo_dict.keys() :
                if key not in edge_dict.keys() :
                    continue
                elif key == 'edge' :
                    attr_d = edge_dict[key]
                    for edge_key,edge_val in obo_dict['edge'].items() :
                        attr_d[edge_key] = edge_val
                else :
                    attr_d = edge_dict[key]
                    val = obo_dict[key]
                    if id not in attr_d.keys() :
                        attr_d[id] = val
        else :
            continue
    ####### save obo information to graph ############
    id = list(edge_dict['name'].keys())
    for go in id :
        name = edge_dict['name'][go]
        root = edge_dict['namespace'][go] 
        if root not in graph_dict.keys() :
            continue 
        else :
            g = graph_dict[root]
            if go not in g.nodes():
                g.add_node(go,namespace = root,name = name)
            else :
                g.nodes[go]['name'] = name
                g.nodes[go]['namespace'] = root
    #save edge infor to graph
    for edge_key,edge_attr in edge_dict['edge'].items() :
        n1,n2 = edge_key
        #makesure both GO terms belong to same GO field
        if edge_dict['namespace'][n1] == edge_dict['namespace'][n2] :
            root = edge_dict['namespace'][n1]
            g = graph_dict[root]
            g.add_edge(n1,n2,relationship = edge_attr)

    return graph_dict,edge_dict

def add_node_level(graph,root,verbose = True):
    '''
    level means minimum path to root node. 
    graph : networkx graph object. GO term as node, with name,namespace attribute
    root : string. GO ID of namespace root node
    '''
    root_dict = {'GO:0008150' : 'biological_process',
                 'GO:0005575' : 'cellular_component',
                 'GO:0003674' : 'molecular_function'}

    node_list = list(graph.nodes())
    node_level = np.repeat(-1,len(node_list))
    if root not in root_dict.keys() :
        print("This root %s have no corresponding namespace" % root)
        return

    for idx,node in enumerate(node_list) :

        namespace = graph.nodes[node]['namespace']
        if namespace == root_dict[root] :
            try :
                l = nx.dijkstra_path_length(graph,node,root)
                node_level[idx] = l
            except :
                if verbose :
                    print("This node : %s can't reach the root !!" % node)
        else :
            if verbose :
                print('Root and namespace not match')
                continue
    level_dict = dict(zip(node_list,node_level))
    nx.set_node_attributes(graph,level_dict,'Level')  # type: ignore


def main() :

    single_go_path = '/home/bruce1996/data/GO/single_obo/'
    go_list = os.listdir(single_go_path)
    graph_d,edge_d = obo2graph(single_go_path)

    namespace2go = {'biological_process' : 'GO:0008150' ,
                    'cellular_component' : 'GO:0005575' ,
                    'molecular_function' :'GO:0003674'}
    for namespace in graph_d.keys() :
        root_go = namespace2go[namespace]
        G = graph_d[namespace]
        add_node_level(G,root_go)
    
    bp = graph_d['biological_process']
    with open("",'wb') as f :
        pickle.dump(bp,f)
    f.close()