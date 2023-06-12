library(survival)
library(survminer)
library(stringr)

condition_list = c('liri_GRN_tumor_only_25/','liri_PIN_tumor_only_25/','liri_GRN_hbv_only_35/','liri_PIN_hbv_only_35/')
summary_output_path = "/home/bruce1996/data/LIHC_anomaly_detection/survival_analysis/5year_coxph_summary/"
for (condition in condition_list){
  survival_path = paste0("/home/bruce1996/data/LIHC_anomaly_detection/survival_analysis/5year_survival_input/",condition)
  file_list = list.files(survival_path)
  output_path = paste0("/home/bruce1996/data/LIHC_anomaly_detection/fig/functional_module/survival_km_plot/",condition)
  for (file in file_list) {
    prefix = substr(file,1,nchar(file)-4)
    print(condition)
    print(prefix)
    res = survival_analysis(paste0(survival_path,file),prefix)
    forest_plot = res[1][[1]]
    km_plot = ggsave_workaround(res[2][[1]])
    model_summary = res[3][[1]]
    write.table(model_summary$coefficients,file = paste0(summary_output_path,str_replace(condition,'/','_'),prefix,"_cox_summary.txt"),sep = '\t',quote = F)
    ggsave(filename = paste0(output_path,prefix,"_forest_plot.png"), plot = forest_plot,width = 8, height = 6 * dim(model_summary$coefficients)[1] / 15, dpi = 300)
    ggsave(filename = paste0(output_path,prefix,"_km_plot.png"), plot = km_plot,width = 8, height = 6, dpi = 300)
  }
}

