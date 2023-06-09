% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/func_ssgsea.R
\name{ssgsea_calculate}
\alias{ssgsea_calculate}
\title{ssgsea_calculate}
\usage{
ssgsea_calculate(exprset = exprset,
signatures = signatures, 
data.type="tpm")
}
\arguments{
\item{exprset}{A expression matrix.}

\item{signatures}{A \code{data.frame} or \code{list} of signature genes.}

\item{data.type}{A type of input expression matrix, "tpm", "fpkm", or "count".}

\item{if_normalize}{Logical. Whether to convert the ssgsea score to \link{0, 1}.}
}
\value{
A data.frame containing the results of ssgsea.
}
\description{
A function to conduct ssgsea.
}
\author{
GM.W
}
