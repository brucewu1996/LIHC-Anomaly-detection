library(stringr)

condition_list = c('lihc_GRN_tumor_only_25/','lihc_PIN_tumor_only_25/','lihc_GRN_hbv_only_35/','lihc_PIN_hbv_only_35/')
summary_output_path = "/home/bruce1996/data/LIHC_anomaly_detection/survival_analysis/5year_coxph_summary/"
for (condition in condition_list){
  survival_path = paste0("/home/bruce1996/data/LIHC_anomaly_detection/survival_analysis/5year_survival_input/",condition)
  file_list = list.files(survival_path)
  output_path = paste0("/home/bruce1996/data/LIHC_anomaly_detection/manuscript/Fig5_survival_of_functional_module/",condition)
  for (file in file_list) {
    prefix = substr(file,1,nchar(file)-4)
    print(condition)
    res = survival_analysis(paste0(survival_path,file),prefix)
    forest_plot = res[1][[1]]
    km_plot = ggsave_workaround(res[2][[1]])
    model_summary = res[3][[1]]
    write.table(model_summary$coefficients,file = paste0(summary_output_path,str_replace(condition,'/','_'),prefix,"_cox_summary.txt"),sep = '\t',quote = F)
    if (dim(model_summary$coefficients)[1] > 20){
      ggsave(filename = paste0(output_path,prefix,"_forest_plot.png"), plot = forest_plot,width = 8, height = 12, dpi = 300)
    }else{
      ggsave(filename = paste0(output_path,prefix,"_forest_plot.png"), plot = forest_plot,width = 8, height = 6, dpi = 300) 
    }
    ggsave(filename = paste0(output_path,prefix,"_km_plot.png"), plot = km_plot,width = 8, height = 6, dpi = 300)
  }
}

