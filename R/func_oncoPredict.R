#' @title run_oncoPredict
#' 
#' @description A function to easily run oncoPredict
#' 
#' @param ExprData A expression matrix.
#' @param batchCorrect The method of batch correction. Choices are "eb" (default) for ComBat, "qn" for quantiles normalization or "none" for no homogenization.
#' @param drugs The candidate drugs selected. Use `mytools::candidate_drugs` to pick some.
#' @param powerTransformPhenotype Should the phenotype be power transformed before we fit the regression model? Default to TRUE, set to FALSE if the phenotype is already known to be highly normal.
#' @param removeLowVaryingGenes What proportion of low varying genes should be removed? 20 percent be default.
#' @param minNumSamples How many training and test samples are required. Print an error if below this threshold.
#' @param printOutput Set to FALSE to supress output.
#' @param removeLowVaringGenesFrom Determine method to remove low varying genes. Options are 'homogenizeData' and 'rawData'.
#' 
#' @return A data.frame containing the predicted drug sensitivity based on exprData. 
#' 
#' @export
#' 
#' @example chemo_res <- run_oncoPredict(ExprData,
#'    batchCorrect = 'eb',
#'    powerTransformPhenotype = FALSE,
#'    removeLowVaryingGenes = 0.2,
#'    minNumSamples = 10,
#'    printOutput = TRUE,
#'    removeLowVaringGenesFrom = 'rawData')
#' 
#' @author GM_W
#' 
run_oncoPredict <- function(ExprData,
                            batchCorrect = 'eb',
                            drugs = NULL,
                            powerTransformPhenotype = TRUE,
                            removeLowVaryingGenes = 0.2,
                            minNumSamples = 10,
                            printOutput = TRUE, 
                            removeLowVaringGenesFrom = 'rawData'){
  if(!file.exists("D:/BaiduSyncdisk/Public_Database/R_functions/oncoPredict/metadata/GDSC2.Rdata"))
  {
    GDSC2_Expr <- readRDS("D:/BaiduSyncdisk/Public_Database/R_functions/oncoPredict/Training Data/GDSC2_Expr (RMA Normalized and Log Transformed).rds")
    GDSC2_Res <- readRDS("D:/BaiduSyncdisk/Public_Database/R_functions/oncoPredict/Training Data/GDSC2_Res.rds")
    # GDSC2_Res <- exp(GDSC2_Res)
    save(GDSC2_Expr,GDSC2_Res,file = "D:/BaiduSyncdisk/Public_Database/R_functions/oncoPredict/metadata/GDSC2.Rdata")
  }else{
    load("D:/BaiduSyncdisk/Public_Database/R_functions/oncoPredict/metadata/GDSC2.Rdata")
  }
  oncoPredict::calcPhenotype(trainingExprData = as.matrix(GDSC2_Expr),
                             trainingPtype = as.matrix(GDSC2_Res),
                             testExprData = as.matrix(ExprData),
                             batchCorrect = batchCorrect,  
                             powerTransformPhenotype = powerTransformPhenotype,
                             removeLowVaryingGenes = removeLowVaryingGenes,
                             minNumSamples = minNumSamples,
                             printOutput = printOutput,
                             removeLowVaringGenesFrom = removeLowVaringGenesFrom)
  chemo.res <- read.csv("./calcPhenotype_Output/DrugPredictions.csv",
                       row.names = 1)
  chemo.res <- exp(chemo.res)
  unlink("calcPhenotype_Output",recursive = TRUE)
  return(chemo.res)
}
