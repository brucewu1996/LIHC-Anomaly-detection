# LIHC-Anomaly-detection
## Introduction 
Hepatocellular carcinoma (HCC) is one of the leading contributors to cancer mortality worldwide and one of the few cancers with a rising mortality rate. Multiple risk factors for HCC have been established including cirrhosis and infection of human immunodeficiency virus (HIV), hepatitis B virus (HBV) and hepatitis C virus (HCV). According to the Global Burden of Disease study on HCC, approximately 42% of global incident cases of HCC between 1990 and 2015 were attributed to chronic HBV infection.
This study addresses the challenge of identifying potential candidates for HBV eradication. Utilizing an *ensemble learning model*, 1,000 linear support vector machine classifiers are trained on TCGA-LIHC dataset. The classifiers assess gene expression profiles, with positive training data from hallmark genes associated with HBV clearance. Negative training data consists of randomly sampled non-hallmark gene expressions. 
This approach aims to uncover genes with patterns resembling those involved in viral clearance, offering insights into potential candidates for further investigation in the context of HBV infection and HCC development.

![image](https://github.com/brucewu1996/LIHC-Anomaly-detection/blob/main/overview.png)

## Data resource  
1. TCGA-LIHC gene expression : [TCGA-LIHC](https://portal.gdc.cancer.gov/projects/TCGA-LIHC).
2. ICGC-LIRI gene expreession : [ICGC-LIRI](https://dcc.icgc.org/projects/LIRI-JP).
3. Hallmark gene set from mSigDB : [WIELAND_UP_BY_HBV_INFECTION geneset](https://www.gsea-msigdb.org/gsea/msigdb/human/geneset/WIELAND_UP_BY_HBV_INFECTION.html).
4. scRNA gene information : [GEO GSE149614](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE149614)
## Tutorial 
1. [Data augmentation](https://github.com/brucewu1996/LIHC-Anomaly-detection/blob/main/scripts/data_augmentation.py) & [ensemble learning](https://github.com/brucewu1996/LIHC-Anomaly-detection/blob/main/scripts/ensemble_learning.py) address by in-house scripts.
2. The ensemble learning model evaluation via [ensemble learning model evaluation](https://github.com/brucewu1996/LIHC-Anomaly-detection/blob/main/scripts/ensemble_model_performance.ipynb).
3. Functional module interaction network construct by [functional module interaction network](https://github.com/brucewu1996/LIHC-Anomaly-detection/blob/main/scripts/functional_module_interaction_network.ipynb).
4. Survival analysis performed by following scripts : [sciprt1](https://github.com/brucewu1996/LIHC-Anomaly-detection/blob/main/scripts/survival_utils.R), [script2](https://github.com/brucewu1996/LIHC-Anomaly-detection/blob/main/scripts/survival_demo.R).
5. Single cell RNA dataset analysis by following : [UMAP clustering](https://github.com/brucewu1996/LIHC-Anomaly-detection/blob/main/scripts/scRNA_candidated_genes_UMAP_clustering.ipynb) & [clustering result evaluation](https://github.com/brucewu1996/LIHC-Anomaly-detection/blob/main/scripts/scRNA_clustering_measurement.ipynb).