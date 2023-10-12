library(stringr)

survival_input_path = "/home/bruce1996/data/LIHC_anomaly_detection/survival_analysis/tcga_lihc_survival_input/"
survival_fig_output_path = "/home/bruce1996/data/LIHC_anomaly_detection/survival_analysis/tcga_lihc_survival_figure/"
survival_cox_output_path = "/home/bruce1996/data/LIHC_anomaly_detection/survival_analysis/tcga_lihc_survival_cox_summary/"
condition_l = c()
ppi_l = c()
module_l = c()
pv_l = c()

for (condition in c('HBV','non-HBV')){
  for (ppi in c('GRN','PIN')){
    survival_path = paste0(survival_input_path,condition,'/',ppi,'/')
    file_list = list.files(survival_path)
    output_path = paste0(survival_fig_output_path,condition,'/',ppi,'/')
    for (file in file_list) {
      prefix = substr(file,1,nchar(file)-4)
      condition_l = append(condition_l,condition)
      ppi_l = append(ppi_l,ppi)
      module_l = append(module_l,prefix)
      res = survival_analysis(paste0(survival_path,file),prefix)
      forest_plot = res[1][[1]]
      logrank_pv = res[3][[1]]$logtest[[3]]
      print(logrank_pv)
      pv_l = append(pv_l,logrank_pv)
      km_plot = ggsave_workaround(res[2][[1]])
      model_summary = res[3][[1]]
      write.table(model_summary$coefficients,file = paste0(survival_cox_output_path,condition,'/',ppi,'/',prefix,"_cox_summary.txt"),sep = '\t',quote = F)
      if (dim(model_summary$coefficients)[1] > 20){
        ggsave(filename = paste0(output_path,prefix,"_forest_plot.png"), plot = forest_plot,width = 8, height = 12, dpi = 300)
      }else{
        ggsave(filename = paste0(output_path,prefix,"_forest_plot.png"), plot = forest_plot,width = 8, height = 6, dpi = 300) 
      }
      ggsave(filename = paste0(output_path,prefix,"_hbv_km_plot.png"), plot = km_plot,width = 8, height = 6, dpi = 300)
      
      df = data.frame(module_l,pv_l,ppi_l,condition_l)
      colnames(df) = c("Module",'pvalue','PPI','HBV')
    }
  }
}

