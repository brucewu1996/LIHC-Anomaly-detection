library(survival)
library(survminer)
data(myeloma)

variable = colnames(myeloma)
cox_variable = variable[! variable %in% c('chr1q21_status', 'event', 'time','molecular_group','treatment')]
formula = as.formula(paste("Surv(time,event) ~",paste(cox_variable, collapse=" + ")))
cox_sex = coxph(formula, data=myeloma)

