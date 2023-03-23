#' @title TIMER_deconv
#'
#' @description Deconvolute bulk expression matrix and predict immune cell proportions.
#' 
#' @usage TIMER_deconv(exprset,
#' indication = "PAAD")
#'
#' @param exprset the expression matrix (not log2 tranformed).
#' @param indication the tumor type (TCGA) of samples.
#' 
#' @return A data.frame containing predicted immune cell proportions.
#' 
#' @export
#' 
#' @example TIMER_deconv(exprset,
#' indication = "PAAD")
#' 
#' @author GM.W
#' 
TIMER_deconv <- function(exprset,indication){
  res_timer <- immunedeconv::deconvolute(exprset,method = "timer",
                                         indications = rep(indication,ncol(exprset)))
  res_timer <- as.data.frame(res_timer)
  rownames(res_timer) <- res_timer$cell_type
  res_timer <- res_timer[,-1]
  res_timer <- as.data.frame(t(res_timer))
  return(res_timer)
}
