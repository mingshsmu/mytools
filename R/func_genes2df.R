#' @title genes2df
#' 
#' @description  A function to convert multi-vectors to df
#' 
#' @param ... The vector of genes.
#' 
#' @return A data.frame
#' 
#' @export
#' 
#' @example genes2df(BER, CPF, FA)
#' colnames(df) = stringr::str_split("BER, CPF, FA",", ")[[1]]
#' CPF = df$CPF[which(df$CPF!=")]
#' 
#' @author GM.W
#' 
genes2df <- function(...){
  list1 = list(...)
  n = length(list1)
  max_genes = max(sapply(list1, length))
  df = data.frame(NA)
  for (i in 1:n) {
    list1[[i]] = c(list1[[i]], rep("",max_genes-length(list1[[i]])))
    df = cbind(df,list1[[i]])
  }
  df = df[,-1,drop=F]
  return(df)
}


