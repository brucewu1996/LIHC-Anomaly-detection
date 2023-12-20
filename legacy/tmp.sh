#!/bin/bash
input_folder="/home/bruce1996/data/LIHC_anomaly_detection/manifold_transformation/differential_expression/full_gene/"
output_folder="/home/bruce1996/data/LIHC_anomaly_detection/manifold_transformation/gsea_summary/"
prefix_list="3A_against_3B 3B_against_4_1 4_1_against_4_2"

for prefix in $prefix_list
do
    echo $prefix
    deseq_path="${input_folder}${prefix}_deseq_result.txt"
    python utils/gsea_for_deseq.py -d $deseq_path -p $prefix -o $output_folder
done