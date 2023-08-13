#' @title TIDE_prepare
#' 
#' @description A function to prepare the input for TIDEpy
#' 
#' @usage TIDE_prepare(exprset,
#' keytype="SYMBOL")
#' 
#' @param exprset genes(row) * samples(column) data.frame
#' @param keytype The type of gene ID. Choose "SYMBOL"(default) or "ENTREZID".
#' 
#' @return A data.frame with ENTREZID rownames
#' 
TIDE_prepare <- function(exprset,keytype="SYMBOL"){

  require(org.Hs.eg.db)
  gene_entrez <- mapIds(x=org.Hs.eg.db, keys = rownames(exprset),
                        keytype = keytype,column = "ENTREZID")
  exprset <- cbind(entrez = gene_entrez,exprset)
  exprset <- na.omit(exprset)
  exprset <- dplyr::distinct(exprset,entrez,.keep_all = T)
  exprset2 <- exprset[,-1]
  rownames(exprset2) <- exprset$entrez
  return(exprset2)
}

#' @title TIDE_pred
#' 
#' @description A function to predict the immunotherapy sensitivity using TIDE algorithm
#' 
#' @usage TIDE_pred(exprset,
#' keytype="SYMBOL")
#' 
#' @param exprset genes(row) * samples(column) data.frame
#' @param keytype The type of gene ID. Choose "SYMBOL"(default) or "ENTREZID".
#' 
#' @return A data.frame.
#' 
#' @export
#' 
#' @author GM.W
#' 
TIDE_pred <- function(exprset, keytype = "SYMBOL"){
  if(!keytype %in% c("SYMBOL", "ENTREZID")){
    stop("Please input the correct keytype! SYMBOL or ENTREZID.")
  }
  
  # ID translation
  if(keytype == "SYMBOL"){
    exprset <- TIDE_prepare(exprset, keytype)
  }
  # write in local for python
  write.csv(exprset, file = "TIDE_res.csv")
  
  # python
  py_code <- c("import pandas as pd", "from tidepy.pred import TIDE")
  py_code <- c(py_code, "df = pd.read_csv('TIDE_res.csv',index_col=0)")
  py_code <- c(py_code, "df = TIDE(df,cancer='Other',pretreat = False,vthres = 0)")
  py_code <- c(py_code, "df.to_csv('TIDE_res.csv')")
  readr::write_lines(x = py_code, file = "py_code.py")
  reticulate::use_condaenv("py3")
  reticulate::source_python(file = "py_code.py")
  
  # Return
  tide <- read.csv("TIDE_res.csv", row.names = 1)
  file.remove("TIDE_res.csv")
  file.remove("py_code.py")
  return(tide)
}

