#' @title fpkm2tpm
#' 
#' @description To convert a fpkm matrix to a tpm matrix
#' 
#' @usage as.data.frame(apply(fpkm, 2, fpkm2tpm))
#' 
#' @param fpkm The expression matrix in the format of `FPKM`.
#' 
#' @return A matrix of TPM.
#' 
#' @export
#' 
#' @examples as.data.frame(apply(fpkm, 2, fpkm2tpm))
#' 
#' @author Not Me.
#' 
fpkm2tpm <- function(fpkm)
{
  exp(log(fpkm) - log(sum(fpkm)) + log(1e6))
}
