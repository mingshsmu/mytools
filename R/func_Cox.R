#' @title Univ_Cox
#' 
#' @description A function to conduct univariate Cox regression
#' 
#' @usage Univ_Cox(factor_list = factor_list, 
#' data_surv = data_surv
#' time = "OS.time",
#' event = "OS.event")
#' 
#' @param factor_list the factor of interset
#' @param data_surv the df containing factor_list & time&event
#' 
#' @return A data.frame contains `p.value`, `HR`, `HR.confint.lower`, `HR.confint.upper`, `HR.CI`.
#' 
#' @export
#' 
#' @example #' @example data_surv <- data.frame(KLF5 = log2(sample(x=1:100, size = 100)+1),
#' PTK7 = log2(sample(x=50:100, size = 100,replace = T)+1),
#' OS.time = sample(x=20:100, size = 100, replace = T),
#' OS.event = sample(x=c(0,1),size = 100, replace = T))
#' factor_list <- c("KLF5","PTK7")
#' Univ_Cox(factor_list = factor_list, data_surv = data_surv)
#' 
#' @author GM.W
#' 
Univ_Cox <- function(factor_list,data_surv,time="OS.time",event="OS.event"){
  require(survival)
  factor_list <-  gsub(pattern = "-",replacement = "..",factor_list)
  factor_list <-  gsub(pattern = "_",replacement = "....",factor_list)
  factor_list <-  gsub(pattern = "/",replacement = "......",factor_list)
  colnames(data_surv) <-  gsub(pattern = "-",replacement = "..",colnames(data_surv))
  colnames(data_surv) <-  gsub(pattern = "_",replacement = "....",colnames(data_surv))
  colnames(data_surv) <-  gsub(pattern = "/",replacement = "......",colnames(data_surv))
  univ_formulas <- sapply(factor_list,
                          function(x) as.formula(paste('Surv(get(time),get(event))~', x)))
  univ_models <- lapply(univ_formulas, function(x){coxph(x, data = data_surv)})
  univ_results <- lapply(univ_models,
                         function(x){
                           x <- summary(x)
                           # p-value
                           p.value<-signif(x$wald["pvalue"], digits=3)
                           # HR
                           HR <-round(x$coef[2], 2);
                           # HR.confidence interval
                           HR.confint.lower <- round(x$conf.int[,"lower .95"], 2)
                           HR.confint.upper <- round(x$conf.int[,"upper .95"],2)
                           HR.CI <- paste0(HR.confint.lower," - ",HR.confint.upper)
                           res<-c(p.value, HR, HR.confint.lower, HR.confint.upper, HR.CI)
                           names(res)<-c("pvalue","Hazard_Ratio","lower","upper","CI")
                           return(res)
                         })
  # result_table
  uni_cox_table <- data.frame(t(as.data.frame(univ_results, check.names = FALSE)))
  rownames(uni_cox_table) = gsub(pattern = "[.]{6}",replacement = "/",rownames(uni_cox_table))
  rownames(uni_cox_table) = gsub(pattern = "[.]{4}",replacement = "_",rownames(uni_cox_table))
  rownames(uni_cox_table) = gsub(pattern = "[.]{2}",replacement = "-",rownames(uni_cox_table))
  uni_cox_table$variable <- rownames(uni_cox_table)
  uni_cox_table$time <- time
  uni_cox_table$event <- event
  uni_cox_table$pvalue = as.numeric(uni_cox_table$pvalue)
  return(uni_cox_table)
}


#' @title Multi_Cox
#' 
#' @description A function to conduct multivariate Cox regression
#' 
#' @usage Multi_Cox(factor_list = factor_list, 
#' data_surv = data_surv,
#' time = "OS.time",
#' event = "OS.event")
#' 
#' @param factor_list the factor of interset
#' @param data_surv the df containing factor_list, time, & event
#' 
#' @return A data.frame contains `p.value`, `HR`, `HR.confint.lower`, `HR.confint.upper`, `HR.CI`.
#' 
#' @export
#' 
#' @example data_surv <- data.frame(KLF5 = log2(sample(x=1:100, size = 100)+1),
#' PTK7 = log2(sample(x=50:100, size = 100,replace = T)+1),
#' OS.time = sample(x=20:100, size = 100, replace = T),
#' OS.event = sample(x=c(0,1),size = 100, replace = T))
#' factor_list <- c("KLF5","PTK7")
#' Multi_Cox(factor_list = factor_list, data_surv = data_surv)
#' 
#' @author GM.W
#' 
Multi_Cox <- function(factor_list, data_surv, time="OS.time",event="OS.event"){
  require(survival)
  factor_list <-  gsub(pattern = "-",replacement = "..",factor_list)
  factor_list <-  gsub(pattern = "_",replacement = "....",factor_list)
  factor_list <-  gsub(pattern = "/",replacement = "......",factor_list)
  colnames(data_surv) <-  gsub(pattern = "-",replacement = "..",colnames(data_surv))
  colnames(data_surv) <-  gsub(pattern = "_",replacement = "....",colnames(data_surv))
  colnames(data_surv) <-  gsub(pattern = "/",replacement = "......",colnames(data_surv))
  myformula1 = as.formula(paste0("Surv(get(time),get(event))~",
                                 paste(factor_list,sep = "",collapse = "+")))
  multi_cox = coxph(myformula1,data = data_surv)
  x = summary(multi_cox)
  # HR
  HR = round(as.matrix(x$coefficients)[,2],2)
  # p-value
  HRp = signif(as.matrix(x$coefficients)[,5],3)
  #HRp = ifelse(HRp < 0.001,"<0.001",round(HRp,3))
  # HR.confidence interval
  HRCILL = round(x$conf.int[,3],2)
  HRCIUL = round(x$conf.int[,4],2)
  CI = paste0(HRCILL," - ",HRCIUL)
  multi_cox_table = as.data.frame(cbind(pvalue=HRp,Hazard_Ratio = HR,lower = HRCILL,upper = HRCIUL,CI=CI))
  rownames(multi_cox_table) = gsub(pattern = "[......]",replacement = "/",rownames(multi_cox_table))
  rownames(multi_cox_table) = gsub(pattern = "[....]",replacement = "_",rownames(multi_cox_table))
  rownames(multi_cox_table) = gsub(pattern = "[..]",replacement = "-",rownames(multi_cox_table))
  multi_cox_table$variable <- rownames(multi_cox_table)
  multi_cox_table$time <- time
  multi_cox_table$event <- event
  return(multi_cox_table)
}
