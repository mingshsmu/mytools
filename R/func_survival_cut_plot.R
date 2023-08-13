#' @title survival_cut_plot
#'
#' @description A function to draw K-M curve according to the median or best cutoff value
#' 
#' @usage survival_cut_plot(surv_data, 
#' variable = "NEIL3", method = "cutpoint",
#'    time = "OS.time",event = "OS",
#'    CI = T, title = "TCGA-PDAC")
#' 
#' @param surv_data The data.frame containing time, event, and genes of interest
#' @param variable Gene or observation of interest
#' @param method Choose one from c("cutpoint","median"), default is "cutpoint".
#' @param time The name of time variable.
#' @param event The name of event variable.
#' @param palette The color schemes for the K-M curves. Default is c("#4dbbd5","#e64b35").
#' @param CI Logical. Whether to draw the Confidence interval. Default is `FALSE`. 
#' @param title The title of the plot. Default is `NULL`.
#' 
#' @return A ggplot object.
#' 
#' @export
#' 
#' @examples survival_cut_plot(surv_data, variable = "NEIL3", method = "cutpoint",
#'    time = "OS.time",event = "OS",CI = T, title = "TCGA-PDAC")
#'    
#' @author GM.W
#' 
survival_cut_plot <- function(surv_data, variable, method=c("cutpoint","median"), 
                              time="time", event="event", palette=c("#4dbbd5","#e64b35"),
                              CI=FALSE, title=NULL) {
  ########### main ##########
  surv_data <- surv_data[,c(time,event,variable)]
  surv_data <- surv_data[!is.na(surv_data[,time]) & !is.na(surv_data[,event]),]
  colnames(surv_data) <- c("time","event","variable")
  if(length(method)>1){method <- method[1]}
  if(method == "cutpoint"){
    cat(crayon::green(">>> Performing the survival analysis on",variable,"according to the **best cutpoint**", 
                     "\n"))
    surv_data <- survminer::surv_cutpoint(surv_data, time = "time", event = "event", variables = "variable")
    surv_data <- survminer::surv_categorize(surv_data)
  }else if(method == "median"){
    cat(crayon::green(">>> Performing the survival analysis on",variable,"according to the **median**", 
                     "\n"))
    surv_data[,"variable"] <- ifelse(surv_data[,"variable"] > median(surv_data[,"variable"]), "high", "low")
  }else{
    stop("The input *method* is wrong!")
  }
  surv_data[,"variable"] <- factor(surv_data[,"variable"],levels = c("low","high"))
  # log-rank test
  surv_diff <- survival::survdiff(as.formula("survival::Surv(time,event)~variable"), data = surv_data)
  pval <- 1 - pchisq(surv_diff$chisq, length(surv_diff$n) - 1)
  HR <- (surv_diff$obs[2]/surv_diff$exp[2])/(surv_diff$obs[1]/surv_diff$exp[1])
  up95 <- exp(log(HR) + qnorm(0.975)*sqrt(1/surv_diff$exp[2]+1/surv_diff$exp[1]))
  low95 <- exp(log(HR) - qnorm(0.975)*sqrt(1/surv_diff$exp[2]+1/surv_diff$exp[1]))
  anno <- paste0("Log-rank\np = ",pval,"\nHazard Ratio = ",HR,"\n95% CI: ",low95,"-",up95)
  # K-M plot
  n <- table(surv_data[,"variable"])
  surv_data[,"variable"] <- ifelse(surv_data[,"variable"]=="low",
                                 paste0("Low (N=",n["low"],")"),
                                 paste0("High (N=",n["high"],")"))
  surv_data[,"variable"] <- factor(surv_data[,"variable"],levels = c(paste0("Low (N=",n["low"],")"),
                                                                 paste0("High (N=",n["high"],")")))
  require(ggkm,warn.conflicts = FALSE, quietly = TRUE)
  p <- ggplot2::ggplot(surv_data,ggplot2::aes(time=time,status=event,color=variable,fill=variable))+
    ggkm::geom_km(size=1.2)+
    ggplot2::annotate("text",x=0,y=0.2,label=anno,fontface = 'italic',hjust=0)+
    ggplot2::labs(x="Follow up time (month)", y="Survival probability")+
    survminer::theme_survminer()+
    ggplot2::theme(legend.position = c(0.7,0.9),plot.title = ggplot2::element_text(hjust = 0.5))+
    ggplot2::scale_color_manual(values = palette)+
    ggplot2::scale_fill_manual(values = palette)+
    ggplot2::labs(color=variable,fill=variable)
  if(CI==TRUE){p <- p+ggkm::geom_kmband()}
  if(!is.null(title)){p <- p+ggplot2::ggtitle(title)}
  return(p)
}


#' @title survival_cut_res
#'
#' @description A function to calculate the p-value of survival analysis.
#' 
#' @usage survival_cut_res(surv_data, 
#'    variable = "NEIL3", 
#'    method = "cutpoint",
#'    time = "OS.time",
#'    event = "OS")
#' 
#' @param surv_data The data.frame containing time, event, and genes of interest
#' @param variable Gene or observation of interest
#' @param method Choose one from c("cutpoint","median"), default is "cutpoint".
#' @param time The name of time variable.
#' @param event The name of event variable.
#' 
#' @return A data.frame containing p-value of survival analysis.
#' 
#' @export
#' 
#' @examples survival_cut_plot(surv_data, variable = "NEIL3", method = "cutpoint",
#'    time = "OS.time",event = "OS")
#'    
#' @author GM.W
#' 
survival_cut_res <- function(surv_data, variable, method=c("cutpoint","median"), 
                              time="time", event="event") {
  ########### main ##########
  surv_data <- surv_data[,c(time,event,variable)]
  colnames(surv_data) <- c("time","event","variable")
  if(length(method)>1){method <- method[1]}
  if(method == "cutpoint"){
    # cat(crayon::green(">>> Performing the survival analysis on",variable,"according to the **best cutpoint**", 
    #                  "\n"))
    surv_data <- survminer::surv_cutpoint(surv_data, time = "time", event = "event", variables = "variable")
    surv_data <- survminer::surv_categorize(surv_data)
  }else if(method == "median"){
    # cat(crayon::green(">>> Performing the survival analysis on",variable,"according to the **median**", 
    #                  "\n"))
    surv_data[,"variable"] <- ifelse(surv_data[,"variable"] > median(surv_data[,"variable"]), "high", "low")
  }else{
    stop("The input *method* is wrong!")
  }
  surv_data[,"variable"] <- factor(surv_data[,"variable"],levels = c("low","high"))
  # log-rank test
  surv_diff <- survival::survdiff(as.formula("survival::Surv(time,event)~variable"), data = surv_data)
  pval <- 1 - pchisq(surv_diff$chisq, length(surv_diff$n) - 1)
  HR <- (surv_diff$obs[2]/surv_diff$exp[2])/(surv_diff$obs[1]/surv_diff$exp[1])
  up95 <- exp(log(HR) + qnorm(0.975)*sqrt(1/surv_diff$exp[2]+1/surv_diff$exp[1]))
  low95 <- exp(log(HR) - qnorm(0.975)*sqrt(1/surv_diff$exp[2]+1/surv_diff$exp[1]))
  anno <- paste0("Log-rank\np = ",pval,"\nHazard Ratio = ",HR,"\n95% CI: ",low95,"-",up95)
  res <- data.frame(variable = variable,method=method, pvalue=pval,hazard_ratio=HR,
                    CI.lower=low95, CI.upper=up95, time=time,event=event)
  return(res)
}