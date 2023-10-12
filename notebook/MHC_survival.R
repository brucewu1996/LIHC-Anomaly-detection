library(survival)
library(survminer)

module = "GO:0019885_antigen-processing-and-presentation-of-endogenous-peptide-antigen-via-MHC-class-I.txt"
survival_path = "/home/bruce1996/data/LIHC_anomaly_detection/survival_analysis/tcga_lihc_survival_input/zscore/HBV/PIN/"

survival_input = read.table(paste0(survival_path,module),header = T,row.names = 1,sep = '\t')
survival_input$Status = survival_input$Status == 'Dead'
survival_input$Stage = factor(survival_input$Stage, levels = c("Stage I","Stage II","Stage III","Stage IV","Not Available"))
# cox regression by gene expression
genes = colnames(survival_input)[8:dim(survival_input)[2]]
# create a formula
fmla <- as.formula(paste("Surv(Survival_day, Status) ~", paste(genes, collapse= "+")))
total_fmla <- as.formula(paste("Surv(Survival_day, Status) ~", paste(colnames(survival_input)[3:dim(survival_input)[2]], collapse= "+")))
model <- coxph( fmla,data = survival_input[c(genes,"Survival_day","Status")] )
risk_score = predict(model, survival_input, type="risk")
survival_input$risk = risk_score > 1
survival_input$risk[survival_input$risk == TRUE] = 'High risk'
survival_input$risk[survival_input$risk == FALSE] = 'Low risk'
output_path = "/home/bruce1996/data/LIHC_anomaly_detection/survival_analysis/tcga_candidate_module/"
write.table(survival_input,file = paste0(output_path,module),sep = '\t',quote = F)


km_col='risk'
km_labels=c("Low risk", "High risk")
km_fmla = as.formula(paste0("Surv(Survival_day, Status) ~ ",km_col))
#fit <- do.call(surfit, args = list(formula = fmla, data = survival_input))
km_survival = data.frame(survival_v2)
fit = survminer::surv_fit(formula = km_fmla,data = survival_input)


kmplot = ggsurvplot(fit, data = survival_input,pval = TRUE,
                    conf.int = TRUE,method = 'survdiff',
                    legend.labs = km_labels,
                    # Add risk table
                    risk.table = TRUE,
                    tables.height = 0.2,
                    palette = c("#98EECC", "#FEA1A1"))
