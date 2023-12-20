import numpy as np
import pandas as pd
import os,pickle,json
import argparse
from sklearn.preprocessing import StandardScaler
from sklearn.linear_model import LogisticRegression
from sklearn.ensemble import RandomForestClassifier
from sklearn.model_selection import train_test_split
from sklearn.metrics import accuracy_score,auc, precision_recall_curve,roc_auc_score

def functional_module_performance(exp_profile,metadata,candidate_go,go_2_gene_dict) :
    '''
    exp_profile : dataframe; must be original un-standardized tcga fpkm file,row is gene, columns is sample 
    metadata : pd.series; index = sample, value = Label (Positive / Negative)
    candidate_go : list;  list of GO terms
    go_2_gene_dict : dict; key = GO term, value = list of GO. Product by functional_module_gsea.py
    '''
    
    y = np.where(metadata == 'Positive',1,0)
    print("Number of %d postive samples in %d samples!" % (sum(y),len(y)))
    metric_matrix = np.zeros([len(candidate_go),7])

    for go_idx,go in enumerate(candidate_go) :
        intersect_gene = set(go_2_gene_dict[go]).intersection(exp_profile.index)
        X = exp_profile.loc[intersect_gene,:].to_numpy().T
        ##add re-standardize section
        scaler = StandardScaler()
        x = scaler.fit_transform(X)
        
        x_train, x_test, y_train, y_test = train_test_split(x,y,train_size=0.3,stratify=y)
        clf = LogisticRegression(random_state=0,max_iter=500).fit(x_train, y_train)
        randomforest = RandomForestClassifier(max_depth=5,n_jobs=32).fit(x_train,y_train)
        y_pred = clf.predict(x_test)
        random_pred = randomforest.predict(x_test)
        lr_precision, lr_recall, _ = precision_recall_curve(y_pred, y_test)
        rf_precision, rf_recall, _ = precision_recall_curve(random_pred, y_test)
        
        metric_matrix[go_idx,0] = len(intersect_gene)
        metric_matrix[go_idx,1] = accuracy_score(y_test,y_pred)
        metric_matrix[go_idx,2] = roc_auc_score(y_test,y_pred)
        metric_matrix[go_idx,3] = auc(lr_recall, lr_precision)
        metric_matrix[go_idx,4] = accuracy_score(y_test,random_pred)
        metric_matrix[go_idx,5] = roc_auc_score(y_test,random_pred)
        metric_matrix[go_idx,6] = auc(rf_recall, rf_precision)

    logistic_df = pd.DataFrame(metric_matrix[:,0:4],index=candidate_go,columns=['Gene number','Accuracy','AUC_ROC','AUC_PR'])
    randomforest_df = pd.DataFrame(metric_matrix[:,(0,4,5,6)],index=candidate_go,columns=['Gene number','Accuracy','AUC_ROC','AUC_PR'])

    return logistic_df,randomforest_df

def convert_go_into_descrption(go_list,go_namespace,go_dict) :
    description = []
    for go in go_list :
        field = go_namespace[go]
        descrip = go_dict[field].nodes[go]['name']
        description.append(descrip)
        
    return description

def main() :
    
    parser = argparse.ArgumentParser()
    parser.add_argument("-e", "--exp_profile",help="path of expression profile")
    parser.add_argument("-m", "--metadata",help="path of HBV label")
    parser.add_argument("-g", "--go",help="path of candidate go term")
    parser.add_argument("-d", "--go_2_gene_dict",help="path of candidate go to gene dict")
    parser.add_argument("-o","--output_path",help = 'path of ensemble model output')
    args = parser.parse_args()
    # load essential file
    exp_profile = pd.read_csv(args.exp_profile,sep='\t',index_col=0)
    label = pd.read_csv(args.metadata,sep='\t',index_col=1)
    metadata = label.loc[list(exp_profile.columns),'HBV']
    
    go_df = pd.read_csv(args.go,sep='\t',index_col=0)
    candidate_go = list(go_df.index)
    with open(args.go_2_gene_dict,'rb') as f :
        go_gene_dict = pickle.load(f)
    f.close()
    ## evaluate candidate performance
    logistic_df,randomforest_df = functional_module_performance(exp_profile,metadata,candidate_go,go_gene_dict)
    ## convert GO accession to description
    #import go graph
    graph_dict = dict()
    graph_path = '/home/bruce1996/data/GO/networkx/'
    go_field = ['biological_process','cellular_component','molecular_function']
    for field in go_field :
        file = graph_path + field + '.go'
        with open(file , 'rb') as f:  # notice the r instead of w
            graph_dict[field] = pickle.load(f)
            f.close()

    with open("/home/bruce1996/data/GO/go2namespace.json",'rb') as f :
        go2namespace = json.load(f)
        f.close() 
        
    description = convert_go_into_descrption(candidate_go,go2namespace,graph_dict)
    
    title = []
    for idx in range(logistic_df.shape[0]) :
        des = description[idx] + ' (n=%d)' % logistic_df['Gene number'][idx]
        title.append(des)
        
    logistic_df['Description'] = title
    randomforest_df['Description'] = title
    
    logistic_df.to_csv(args.output_path + 'logistic_performance.txt',sep='\t')
    randomforest_df.to_csv(args.output_path + 'randomforest_performance.txt',sep='\t')
    
if __name__ == '__main__' :
    main()