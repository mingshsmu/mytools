% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/func_TIDE_preparing.R
\name{TIDE_pred}
\alias{TIDE_pred}
\title{TIDE_pred}
\usage{
TIDE_pred(exprset,
keytype="SYMBOL")
}
\arguments{
\item{exprset}{genes(row) * samples(column) data.frame}

\item{keytype}{The type of gene ID. Choose "SYMBOL"(default) or "ENTREZID".}
}
\value{
A data.frame.
}
\description{
A function to predict the immunotherapy sensitivity using TIDE algorithm
}
\author{
GM.W
}
