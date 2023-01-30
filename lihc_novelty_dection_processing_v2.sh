#!/bin/bash
#this script is for whole lihc novelty dection processing

lihc_path="/home/bruce1996/data/LIHC_anomaly_detection"
biomart_path="/home/bruce1996/data/LIHC_anomaly_detection/data/coding_gene_info/biomart_protein_coding_gene.txt"
go_path="/home/bruce1996/data/GO/networkx/"
prefix_list=('tumor_only' 'normal_only' 'hbv_only' 'nonhbv_only')
hallmark="${lihc_path}/data/hallmark_gene/hallmark_protein_coding_ens_id_without_outlier.txt"
for prefix in "${prefix_list[@]}"
do
    project_prefix="${prefix}_std_by_gene_without_outlier"
    echo "LIHC novelty dectection for ${project_prefix}"
    #first step is data augmentation
    exp_profile="${lihc_path}/ensemble_training/${prefix}_coding_gene_${condition}_with_synthetic.txt"
    python data_augmentation.py -i "${lihc_path}/data/exp_profile/lihc_coding_gene_std_by_gene_${prefix}.txt" \
		-r $hallmark \
		-f "/home/bruce1996/data/LIHC_anomaly_detection/fig/data_augmentation/${project_prefix}_" \
		-o $exp_profile \
		-n 176
    log_ts=`date "+%Y-%m-%d %H:%M:%S"`
    echo "Synthetic data in : ${exp_profile}" >> lihc_novelty_associated_files.txt
    echo "$log_ts Data augmentation stage is completed !!" >> lihc_novelty_log.txt
    #second stage is ensemble learning
    echo "Ensemble learning for $prefix section"
    ensemble_result="${lihc_path}/ensemble_result/${project_prefix}/"
	python ensemble_learning.py -e $exp_profile \
			-c $hallmark \
			-r 1000 -t 64 \
			-o $ensemble_result \
			-p $prefix --stand False
    log_ts=`date "+%Y-%m-%d %H:%M:%S"`
    echo "Ensemble learning result in : ${ensemble_result}" >> lihc_novelty_associated_files.txt
    echo "$log_ts Ensemble learning stage is completed !!" >> lihc_novelty_log.txt
    #Step 2.1 ensemble learning result visualization
    ensemble_performance="${lihc_path}/fig/ensemble_model_performance/"
    ensemble_composition="${lihc_path}/fig/candidate_gene_composition/"
    python ensemble_model_performance.py -v $ensemble_result \
                        -p $prefix \
                        -o $ensemble_performance \
                        -f "${project_prefix}_"

    python candidate_gene_characteristic.py -i $exp_profile \
                        -r $hallmark \
                        -v $ensemble_result \
                        -p $prefix \
                        -n 30 \
                        -o $ensemble_composition \
                        -f "${project_prefix}_"
    echo "Ensemble learning performance in : ${ensemble_performance}" >> lihc_novelty_associated_files.txt
    echo "Ensemble learning composition in : ${ensemble_composition}" >> lihc_novelty_associated_files.txt
    echo 'Ensemble learning visualization stage is completed !!' >> lihc_novelty_log.txt
    #Step 3 Fisher exact test
    echo "Fisher exact test for prefix : ${prefix}"
    fisher_result="${lihc_path}/functional_profiling/fisher_exact_test/"
    vote="${lihc_path}/ensemble_result/vote_result/${prefix}_std_by_gene_without_outlier_vote_np_ratio_25.txt"
    python validated_edge_v2.py -i $exp_profile -b $biomart_path -g $go_path -v $vote \
                    -c 0.4 -t 32 \
                    -o $fisher_result -n $project_prefix
    log_ts=`date "+%Y-%m-%d %H:%M:%S"`
    echo "Fisher exact test sparse matrix in : ${fisher_result}${project_prefix}_validated_gene2go_matrix.npz" >> lihc_novelty_associated_files.txt
    echo "Fisher exact test go 2 gene dict in : ${fisher_result}${project_prefix}_validated_go_2_gene_dict.pkl" >> lihc_novelty_associated_files.txt
    echo "$log_ts Fisher exact test stage is completed !!" >> lihc_novelty_log.txt
    #Step 4 Gene enrichment analysis
    echo "Gene enrichment analysis for prefix : ${prefix}"
    gsea_result="${lihc_path}/functional_profiling/candidate_functional_module/"
    python functional_module_gsea.py -g "${fisher_result}${project_prefix}_validated_go_2_gene_dict.pkl" \
                                    -e $exp_profile \
                                    -v $vote \
                                    -d "${lihc_path}/differential_expression/lihc_${prefix}_deseq_result.txt" \
                                    -o $gsea_result \
                                    -p $project_prefix -w True
    log_ts=`date "+%Y-%m-%d %H:%M:%S"`
    echo "GSEA activated functional module in : ${gsea_result}${project_prefix}_activated_functional_module.txt" >> lihc_novelty_associated_files.txt
    echo "GSEA inactivated functional module in : ${gsea_result}${project_prefix}_inactivated_functional_module.txt" >> lihc_novelty_associated_files.txt
    echo "$log_ts Gene set enrichment analysis stage is completed !!" >> lihc_novelty_log.txt
    #Step 5 evaluate model performance 
    echo "Model performance of prefix : ${prefix}"
    hbv_prediction_result="${lihc_path}/functional_profiling/functional_module_evaluation/${project_prefix}/"
    python evaluate_functional_module_performance.py -e "/home/bruce1996/data/LIHC_anomaly_detection/data/exp_profile/LIHC_tumor_coding_gene_fpkm.txt" \
                    -m '/home/bruce1996/data/LIHC_anomaly_detection/data/sample_info/sample_info_df.txt' \
                    -g "${gsea_result}${project_prefix}_activated_functional_module.txt" \
                    -d "${fisher_result}${project_prefix}_validated_go_2_gene_dict.pkl" \
                    -o "${hbv_prediction_result}hbv_prediction_performance_activated_${prefix}_"

    python evaluate_functional_module_performance.py -e "/home/bruce1996/data/LIHC_anomaly_detection/data/exp_profile/LIHC_tumor_coding_gene_fpkm.txt" \
                    -m '/home/bruce1996/data/LIHC_anomaly_detection/data/sample_info/sample_info_df.txt' \
                    -g "${gsea_result}${project_prefix}_inactivated_functional_module.txt" \
                    -d "${fisher_result}${project_prefix}_validated_go_2_gene_dict.pkl" \
                    -o "${hbv_prediction_result}hbv_prediction_performance_inactivated_${prefix}_"

    log_ts=`date "+%Y-%m-%d %H:%M:%S"`
    echo "HBV prediction by activated functional module result in : ${hbv_prediction_result}hbv_prediction_performance_activated_${prefix}_logistic/randomforest_performance.txt" >> lihc_novelty_associated_files.txt
    echo "HBV prediction by inactivated functional module result in : ${hbv_prediction_result}hbv_prediction_performance_inactivated_${prefix}_logistic/randomforest_performance.txt" >> lihc_novelty_associated_files.txt
    echo "$log_ts HBV prediction stage is completed !!" >> lihc_novelty_log.txt
    #step 5.1 plot model performance
    echo "HBV prediction performance for prefix : ${prefix}"
    hbv_prediction_performance="${lihc_path}/fig/functional_module/hbv_prediction_performance/${project_prefix}/"
    python candidate_functional_module.py -v $vote \
                    -g "${gsea_result}${project_prefix}_activated_functional_module.txt" \
                    -d "${fisher_result}${project_prefix}_validated_go_2_gene_dict.pkl" \
                    -l "${hbv_prediction_result}hbv_prediction_performance_activated_${prefix}_logistic_performance.txt" \
                    -r "${hbv_prediction_result}hbv_prediction_performance_activated_${prefix}_randomforest_performance.txt" \
                    -p "activate_${project_prefix}" \
                    -o $hbv_prediction_performance

    python candidate_functional_module.py -v $vote \
                    -g "${gsea_result}${project_prefix}_inactivated_functional_module.txt" \
                    -d "${fisher_result}${project_prefix}_validated_go_2_gene_dict.pkl" \
                    -l "${hbv_prediction_result}hbv_prediction_performance_inactivated_${prefix}_logistic_performance.txt" \
                    -r "${hbv_prediction_result}hbv_prediction_performance_inactivated_${prefix}_randomforest_performance.txt" \
                    -p "inactivate_${project_prefix}" \
                    -o $hbv_prediction_performance

    log_ts=`date "+%Y-%m-%d %H:%M:%S"`
    echo "HBV prediction by activated functional module result visualization in : ${hbv_prediction_performance}hbv_prediction_by_logistic/randomforest_activate_${prefix}.png" >> lihc_novelty_associated_files.txt
    echo "HBV prediction by inactivated functional module result in : ${hbv_prediction_performance}hbv_prediction_by_logistic/randomforest_inactivate_${prefix}.png" >> lihc_novelty_associated_files.txt
    echo "$log_ts HBV prediction visualization stage is completed !!" >> lihc_novelty_log.txt
    
done