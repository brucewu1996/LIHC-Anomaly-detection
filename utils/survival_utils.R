library(survival)
library(survminer)
library(stringr)


survival_analysis<-function(survival_info,fig_prefix,km_col='risk',km_labels=c("Low risk", "High risk")){
  
  survival_input = read.table(survival_info,header = T,row.names = 1,sep = '\t')
  # change status info from str into boolean
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
  survival_v2 = data.frame(survival_input)
  # final cox regression
  model_v2 <- coxph(total_fmla,data = survival_v2)
  forest_plot = ggforest(model_v2,main=fig_prefix,data=survival_v2)
  ## km plot
  km_fmla = as.formula(paste0("Surv(Survival_day, Status) ~ ",km_col))
  #fit <- do.call(surfit, args = list(formula = fmla, data = survival_input))
  km_survival = data.frame(survival_v2)
  fit = survminer::surv_fit(formula = km_fmla,data = km_survival)
  logrank_pv = surv_pvalue(fit)
  kmplot = ggsurvplot(fit, data = km_survival,pval = TRUE,
                      conf.int = TRUE,method = 'survdiff',
                      legend.title = str_to_title(km_col),
                      legend.labs = km_labels,
                      # Add risk table
                      risk.table = TRUE,
                      tables.height = 0.2,
                      palette = c("#98EECC", "#FEA1A1")) + ggtitle(fig_prefix)
  return(list(forest_plot,kmplot,summary(model),logrank_pv))
}

survival_analysis_v2<-function(clinical_info,exp_m,fig_prefix,risk_score_mode=TRUE,confounding_factor = c("Gender","Age","Stage"),
                               km_col='risk',km_labels=c("Low risk", "High risk")){
  
  #clinical_info = read.table(clinical_path,header = T,row.names = 1,sep = '\t')
  clinical_info = clinical_info[,c("Survival_day","Status",confounding_factor)]
  #exp_m = read.table(exp_path,header = T,row.names = 1,sep = '\t')
  # change status info from str into boolean
  clinical_info$Status = survival_input$Status == 'Dead'
  #clinical_info$Stage = factor(survival_input$Stage, levels = c("Stage I","Stage II","Stage III","Stage IV","Not Available"))
  # create a formula for cox regression by gene expression
  genes = colnames(exp_m)
  fmla <- as.formula(paste("Surv(Survival_day, Status) ~", paste(genes,collapse= "+")))
  #merge clinical info & expression matrix
  survival_input = cbind.data.frame(clinical_info,exp_m)
  model <- coxph(fmla,data = survival_input[c(genes,"Survival_day","Status")] )
  risk_score = predict(model, survival_input, type="risk")
  survival_input$risk = risk_score > 1
  survival_input$risk_score = risk_score
  survival_v2 = data.frame(survival_input)
  # Second stage cox regression
  if (risk_score_mode){
    total_fmla <- as.formula(paste("Surv(Survival_day, Status) ~", paste(c(confounding_factor,'risk_score'),collapse= "+")))
  }else {
    total_fmla <- as.formula(paste("Surv(Survival_day, Status) ~", paste(c(confounding_factor,genes),collapse= "+")))
  }
  model_v2 <- coxph(total_fmla,data = survival_v2)
  ## Visualization cox result
  forest_plot = ggforest(model_v2,main=fig_prefix,data=survival_v2)
  km_fmla = as.formula(paste0("Surv(Survival_day, Status) ~ ",km_col))
  #fit <- do.call(surfit, args = list(formula = fmla, data = survival_input))
  km_survival = data.frame(survival_v2)
  fit = survminer::surv_fit(formula = km_fmla,data = km_survival)
  kmplot = ggsurvplot(fit, data = km_survival,pval = TRUE,
                      conf.int = TRUE,
                      legend.title = str_to_title(km_col),
                      legend.labs = km_labels,
                      # Add risk table
                      risk.table = TRUE,
                      tables.height = 0.2,
                      palette = c("#98EECC", "#FEA1A1")) + ggtitle(fig_prefix)
  return(list(forest_plot,kmplot,summary(model)))
}

ggsave_workaround <- function(g){survminer:::.build_ggsurvplot(x = g,
                                                               surv.plot.height = NULL,
                                                               risk.table.height = NULL,
                                                               ncensor.plot.height = NULL)}