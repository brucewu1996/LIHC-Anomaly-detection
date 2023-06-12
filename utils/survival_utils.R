library(survival)
library(survminer)

survival_analysis<-function(survival_info,fig_prefix){
  
  survival_input = read.table(survival_info,header = T,row.names = 1,sep = '\t')
  # change status info from str into boolean
  survival_input$Status = survival_input$Status == 'True'
  survival_input$Stage = factor(survival_input$Stage, levels = c("Stage I","Stage II","Stage III","Stage IV","Not Available"))
  # cox regression by gene expression
  genes = colnames(survival_input)[6:dim(survival_input)[2]]
  # create a formula
  fmla <- as.formula(paste("Surv(Survival_day, Status) ~", paste(genes, collapse= "+")))
  total_fmla <- as.formula(paste("Surv(Survival_day, Status) ~", paste(colnames(survival_input)[3:dim(survival_input)[2]], collapse= "+")))
  model <- coxph( fmla,data = survival_input[c(genes,"Survival_day","Status")] )
  risk_score = predict(model, survival_input, type="risk")
  survival_input$risk = risk_score > 1
  # final cox regression
  model <- coxph(total_fmla,data = survival_input)
  forest_plot = ggforest(model,main = fig_prefix)
  ## km plot
  fmla = as.formula("Surv(Survival_day, Status) ~ risk")
  fit <- do.call(survfit, args = list(formula = fmla, data = survival_input))
  kmplot = ggsurvplot(fit, data = survival_input,pval = TRUE,
                      conf.int = TRUE,
                      legend.title = "Risk",
                      legend.labs = c("Low risk", "High risk"),
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