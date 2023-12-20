#!/bin/bash
lihc_path="/home/bruce1996/data/LIHC_anomaly_detection"
prefix_list='tumor_only normal_only hbv_only nonhbv_only'   
for prefix in $prefix_list
do
    echo "prefix is ${prefix}"

    python ensemble_model_performance.py -v "/home/bruce1996/data/LIHC_anomaly_detection/ensemble_result/${prefix}_std_by_sample/" \
                        -p "${prefix}" \
                        -o "/home/bruce1996/data/LIHC_anomaly_detection/fig/ensemble_model_performance/" \
                        -f "${prefix}_std_by_sample"

    python candidate_gene_characteristic.py -i "${lihc_path}/ensemble_training/${prefix}_coding_gene_std_by_sample_with_synthetic.txt" \
                        -r "${lihc_path}/data/hallmark_gene/hallmark_protein_coding_ens_id.txt" \
                        -v "/home/bruce1996/data/LIHC_anomaly_detection/ensemble_result/${prefix}_std_by_sample/" \
                        -p "${prefix}" \
                        -n 25 \
                        -o "/home/bruce1996/data/LIHC_anomaly_detection/fig/candidate_gene_composition/" \
                        -f "${prefix}_std_by_sample"
done