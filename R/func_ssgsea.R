#' @title ssgsea_calculate
#'
#' @description  A function to conduct ssgsea.
#' 
#' @usage ssgsea_calculate(exprset = exprset,
#' signatures = signatures, 
#' data.type="tpm")
#' 
#' @param exprset A expression matrix.
#' @param signatures A `data.frame` or `list` of signature genes.
#' @param data.type A type of input expression matrix, "tpm", "fpkm", or "count".
#' @param if_normalize Logical. Whether to convert the ssgsea score to [0, 1].
#' 
#' @return A data.frame containing the results of ssgsea.
#' 
#' @export
#' 
#' @example ssgsea_calculate(exprset = exprset,
#' signatures = signatures, 
#' data.type="tpm")
#' 
#' @author GM.W
#' 
ssgsea_calculate <- function(exprset, signatures, 
                             data.type=c("tpm","fpkm","count"),
                             if_normalize = FALSE){
  
  #' @title .normalize
  .normalize <- function(x){
    return((x - min(x)) / (max(x)-min(x)))
  }
  
  message(">>> Initializing...")
  data.type <- data.type[1]
  
  # fpkm or tpm
  if(data.type %in% c("fpkm","tpm")){
    if(max(exprset) > 100){
      # log2(x+1) transform if necessary
      message(">>> The tpm/fpkm expression matrix has not been log2(x+1) tranformed.")
      message(">>> Performing log2(x+1) tranforming...")
      exprset <- log2(exprset+1)
    }
    # ssgsea
    res <- GSVA::gsva(as.matrix(exprset),signatures,method='ssgsea',
                      kcdf='Gaussian',abs.ranking=TRUE,verbose=T)
  }
  
  # count
  if(data.type == "count"){
    # round(2^x-1) if necessary 
    if(max(exprset) < 100){
      message(">>> The count expression matrix has been log2(x+1) tranformed.")
      message(">>> Performing round(2^x-1) tranforming...")
      exprset <- round(2^exprset-1)
    }
    res <- GSVA::gsva(as.matrix(exprset),signatures,method='ssgsea',
               kcdf='Poisson',ssgsea.norm=T,abs.ranking=TRUE,verbose=T)
  }
  res <- t(res)
  
  # if_normalize
  if(if_normalize){
    res <- apply(res, 2, .normalize)
  }
  res <- as.data.frame(res)
  message("***Done!***")
  return(res)
}