#!/bin/bash
<<comment
this script is for whole lihc novelty dection processing
comment

lihc_path="/home/bruce1996/data/LIHC_anomaly_detection"
prefix_list='tumor_only normal_only hbv_only nonhbv_only' 
#first step is data augmentation
<<comment
for prefix in $prefix_list
do
	python data_augmentation.py -i "${lihc_path}/lihc_coding_gene_std_by_sample_${prefix}.txt" \
		-r "/home/bruce1996/data/LIHC_anomaly_detection/data/hallmark_gene/hallmark_protein_coding_ens_id.txt" \
		-f "/home/bruce1996/data/LIHC_anomaly_detection/fig/data_augmentation/${prefix}_std_by_sample_" \
		-o "/home/bruce1996/data/LIHC_anomaly_detection/ensemble_training/${prefix}_coding_gene_std_by_sample_with_synthetic.txt" \
		-n 198
done
echo 'data augmentation stage is completed !!' >> log.txt

#second is ensemble learning
for prefix in $prefix_list
do
	echo "Ensemble learning for $prefix section"
	python ensemble_learning.py -e "${lihc_path}/data/exp_profile/lihc_coding_gene_std_by_sample_${prefix}.txt" \
			-c "${lihc_path}/data/hallmark_gene/hallmark_protein_coding_ens_id.txt" \
			-r 1000 -t 64 \
			-o "${lihc_path}/ensemble_result/${prefix}_std_by_sample/" \
			-p $prefix --stand False
done
echo 'Ensemble learning stage is completed !!' >> log.txt
#Step 2.1 ensemble learning result visualization
for prefix in $prefix_list
do
    echo "prefix is ${prefix}"

    python ensemble_model_performance.py -v "/home/bruce1996/data/LIHC_anomaly_detection/ensemble_result/${prefix}_std_by_sample/" \
                        -p "${prefix}" \
                        -o "/home/bruce1996/data/LIHC_anomaly_detection/fig/ensemble_model_performance/" \
                        -f "${prefix}_std_by_sample_"

    python candidate_gene_characteristic.py -i "${lihc_path}/data/exp_profile/lihc_coding_gene_std_by_sample_${prefix}.txt" \
                        -r "${lihc_path}/data/hallmark_gene/hallmark_protein_coding_ens_id.txt" \
                        -v "/home/bruce1996/data/LIHC_anomaly_detection/ensemble_result/${prefix}_std_by_sample/" \
                        -p "${prefix}" \
                        -n 25 \
                        -o "/home/bruce1996/data/LIHC_anomaly_detection/fig/candidate_gene_composition/" \
                        -f "${prefix}_std_by_sample_"
done

#Step 3 Fisher exact test
for prefix in $prefix_list
do
    echo "prefix is ${prefix}"

    python validated_edge_v2.py -i "${lihc_path}/data/exp_profile/lihc_coding_gene_std_by_sample_${prefix}.txt" \
                    -b "/home/bruce1996/data/LIHC_anomaly_detection/data/coding_gene_info/biomart_protein_coding_gene.txt" \
                    -g "/home/bruce1996/data/GO/networkx/" \
                    -v "${lihc_path}/ensemble_result/vote_result/lihc_ensemble_vote_result_${prefix}_std_by_sample_np_25.txt" \
                    -c 0.4 -t 32 \
                    -o "${lihc_path}/functional_profiling/fisher_exact_test/" \
                    -n "${prefix}_std_by_sample"

done
echo 'Fisher exact test stage is completed !!' >> log.txt
#Step 4 GSEA
for prefix in $prefix_list
do
    echo "prefix is ${prefix}"
    python functional_module_gsea.py -g "${lihc_path}/functional_profiling/fisher_exact_test/${prefix}_std_by_sample_validated_go_2_gene_dict.pkl" \
                                    -v "${lihc_path}/ensemble_result/vote_result/lihc_ensemble_vote_result_${prefix}_std_by_sample_np_25.txt" \
                                    -d "${lihc_path}/differential_expression/lihc_${prefix}_deseq_result.txt" \
                                    -o "${lihc_path}/functional_profiling/candidate_functional_module/" \
                                    -p "${prefix}_std_by_sample"
done
echo 'GSEA stage is completed !!' >> log.txt


#Step 5 evaluate model performance 
for prefix in $prefix_list
do
    echo "prefix is ${prefix}"

    python evaluate_functional_module_performance.py -e "/home/bruce1996/data/LIHC_anomaly_detection/data/exp_profile/lihc_protein_coding_gene_std_exp_profile_tumor_only.txt" \
                    -m '/home/bruce1996/data/LIHC_anomaly_detection/data/sample_info/sample_info_df.txt' \
                    -g "/home/bruce1996/data/LIHC_anomaly_detection/functional_profiling/candidate_functional_module/${prefix}_std_by_sample_inactivated_functional_module.txt" \
                    -d "/home/bruce1996/data/LIHC_anomaly_detection/functional_profiling/fisher_exact_test/${prefix}_std_by_sample_validated_go_2_gene_dict.pkl" \
                    -o "/home/bruce1996/data/LIHC_anomaly_detection/functional_profiling/functional_module_evaluation/hbv_prediction_performance_inactivated_${prefix}_"
    
    python evaluate_functional_module_performance.py -e "/home/bruce1996/data/LIHC_anomaly_detection/data/exp_profile/lihc_protein_coding_gene_std_exp_profile_tumor_only.txt" \
                    -m '/home/bruce1996/data/LIHC_anomaly_detection/data/sample_info/sample_info_df.txt' \
                    -g "/home/bruce1996/data/LIHC_anomaly_detection/functional_profiling/candidate_functional_module/${prefix}_std_by_sample_activated_functional_module.txt" \
                    -d "/home/bruce1996/data/LIHC_anomaly_detection/functional_profiling/fisher_exact_test/${prefix}_std_by_sample_validated_go_2_gene_dict.pkl" \
                    -o "/home/bruce1996/data/LIHC_anomaly_detection/functional_profiling/functional_module_evaluation/hbv_prediction_performance_activated_${prefix}_"

done
comment
#step 5.1 plot model performance
for prefix in $prefix_list
do
    echo "prefix is ${prefix}"
    #activate functional module
    python candidate_functional_module.py -v "${lihc_path}/ensemble_result/vote_result/lihc_ensemble_vote_result_${prefix}_std_by_sample_np_25.txt" \
                    -g "/home/bruce1996/data/LIHC_anomaly_detection/functional_profiling/candidate_functional_module/${prefix}_std_by_sample_activated_functional_module.txt" \
                    -d "/home/bruce1996/data/LIHC_anomaly_detection/functional_profiling/fisher_exact_test/${prefix}_std_by_sample_validated_go_2_gene_dict.pkl" \
                    -l "${lihc_path}/functional_profiling/functional_module_evaluation/hbv_prediction_performance_activated_${prefix}_logistic_performance.txt" \
                    -r "${lihc_path}/functional_profiling/functional_module_evaluation/hbv_prediction_performance_activated_${prefix}_randomforest_performance.txt" \
                    -p "activate_${prefix}" \
                    -o "${lihc_path}/fig/functional_module/hbv_prediction_performance/"
    #inactivate functional module
    python candidate_functional_module.py -v "${lihc_path}/ensemble_result/vote_result/lihc_ensemble_vote_result_${prefix}_std_by_sample_np_25.txt" \
                    -g "/home/bruce1996/data/LIHC_anomaly_detection/functional_profiling/candidate_functional_module/${prefix}_std_by_sample_inactivated_functional_module.txt" \
                    -d "/home/bruce1996/data/LIHC_anomaly_detection/functional_profiling/fisher_exact_test/${prefix}_std_by_sample_validated_go_2_gene_dict.pkl" \
                    -l "${lihc_path}/functional_profiling/functional_module_evaluation/hbv_prediction_performance_inactivated_${prefix}_logistic_performance.txt" \
                    -r "${lihc_path}/functional_profiling/functional_module_evaluation/hbv_prediction_performance_inactivated_${prefix}_randomforest_performance.txt" \
                    -p "inactivate_${prefix}" \
                    -o "${lihc_path}/fig/functional_module/hbv_prediction_performance/"
done
