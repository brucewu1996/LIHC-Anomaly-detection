library(survival)
library(survminer)

cox_input_file = '/home/bruce1996/data/LIHC_anomaly_detection/manuscript/material/module_gene_hbv_prediction/module_single_sample_gsea_zscore_cox_regression_input.txt'
survival_input = read.table(cox_input_file,header = T,row.names = 1,sep = '\t')
survival_input$Status = survival_input$Status == 'True'
survival_input$Stage = factor(survival_input$Stage, levels = c("Stage I","Stage II","Stage III","Stage IV","Not Available"))
# cox regression by gene expression
features = colnames(survival_input)[7:dim(survival_input)[2]]
total_fmla <- as.formula(paste("Surv(Survival_day, Status) ~", paste(colnames(survival_input)[3:dim(survival_input)[2]], collapse= "+")))
model_v2 <- coxph(total_fmla,data = survival_input)
forest_plot = ggforest(model_v2,main='Single sample GSEA',data=survival_input)
output_path = '/home/bruce1996/data/LIHC_anomaly_detection/manuscript/material/module_gene_hbv_prediction/'
ggsave(filename = paste0(output_path,'single_sample_gsea',"_forest_plot.png"), plot = forest_plot,width = 12, height =15 , dpi = 300)