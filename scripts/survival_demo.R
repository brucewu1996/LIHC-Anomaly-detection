# Load requirement
survival_input = read.table("/home/bruce1996/repo/LIHC_anomaly_detection/data/Survival-analysis-exp-matrix/lihc_GRN_hbv_only_35/GO:0045824_negative-regulation-of-innate-immune-response.txt",header = T,row.names = 1,sep = '\t')
confounding_factor = c("Age")
gene_list = colnames(survival_input)[8:dim(survival_input)[2]]
exp_m = survival_input[,gene_list]
clinical_info = survival_input[,c("Survival_day","Status","Gender","Age","Stage")]
# Parameter setting
fig_title="negative-regulation-of-innate-immune-response"
prefix="GO:0045824_negative-regulation-of-innate-immune-response"
survival_fig_output_path="/home/bruce1996/repo/LIHC_anomaly_detection/data/test_survival/"
survival_cox_output_path="/home/bruce1996/repo/LIHC_anomaly_detection/data/test_survival/"
res = survival_analysis_v2(clinical_info,exp_m,fig_title,
                           confounding_factor = c("Age","Gender"),
                           risk_score_mode = FALSE)
#Save the survival result into different variable
forest_plot = res[1][[1]]
km_plot = ggsave_workaround(res[2][[1]])
model_summary = round(res[3][[1]]$coefficients,4)
## output the result
write.table(model_summary,file = paste0(survival_cox_output_path,prefix,"_cox_summary.txt"),sep = '\t',quote = F)
if (dim(model_summary)[1] > 20){
  ggsave(filename = paste0(survival_fig_output_path,prefix,"_forest_plot.pdf"), plot = forest_plot,width = 10, height = 12, dpi = 300,device = 'pdf')
}else{
  ggsave(filename = paste0(survival_fig_output_path,prefix,"_forest_plot.pdf"), plot = forest_plot,width = 10, height = 6, dpi = 300,device = 'pdf') 
}
ggsave(filename = paste0(survival_fig_output_path,prefix,"_hbv_km_plot.pdf"), plot = km_plot,width = 8, height = 6, dpi = 300,device = 'pdf')

