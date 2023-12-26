library(stringr)

input_path = "/mnt/sdb/bruce/LIHC_anomaly_detection/survival_analysis/survival_exp_matrix/zscore/"
fig_output_path = "/mnt/sdb/bruce/LIHC_anomaly_detection/survival_analysis/survival_figure/zscore/"
cox_output_path = "/mnt/sdb/bruce/LIHC_anomaly_detection/survival_analysis/survival_cox_summary/zscore/"
for (folder in dir(input_path)){
  survival_path = paste0(input_path,folder,'/')
  survival_fig_output_path = paste0(fig_output_path,folder,'_risk_score/')
  survival_cox_output_path = paste0(cox_output_path,folder,'_risk_score/')
  # if output folder no exist, create folder
  if (dir.exists(survival_fig_output_path) == FALSE){
    dir.create(survival_fig_output_path)
  }
  if (dir.exists(survival_cox_output_path) == FALSE){
    dir.create(survival_cox_output_path)
  }
  for (module in list.files(survival_path)){
    print(paste0("survival analysis for module :",survival_path,module))
    survival_input = read.table(paste0(survival_path,module),header = T,row.names = 1,sep = '\t')
    gene_list = colnames(survival_input)[8:dim(survival_input)[2]]
    confounding_factor = c("Gender","Age","Stage")
    clinical_info = survival_input[,c("Survival_day","Status",confounding_factor)]
    exp_m = survival_input[,gene_list]
    prefix = substr(module,1,nchar(module)-4)
    res = survival_analysis_v2(clinical_info,exp_m,strsplit(prefix,'_')[[1]][2])
    forest_plot = res[1][[1]]
    km_plot = ggsave_workaround(res[2][[1]])
    model_summary = round(res[3][[1]]$coefficients,4)
    
    write.table(model_summary,file = paste0(survival_cox_output_path,prefix,"_cox_summary.txt"),sep = '\t',quote = F)
    if (dim(model_summary)[1] > 20){
      ggsave(filename = paste0(survival_fig_output_path,prefix,"_forest_plot.png"), plot = forest_plot,width = 10, height = 12, dpi = 300)
      }else{
      ggsave(filename = paste0(survival_fig_output_path,prefix,"_forest_plot.png"), plot = forest_plot,width = 10, height = 6, dpi = 300) 
    }
    ggsave(filename = paste0(survival_fig_output_path,prefix,"_hbv_km_plot.png"), plot = km_plot,width = 8, height = 6, dpi = 300)
  }
}

