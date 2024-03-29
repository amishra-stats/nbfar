% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/misc.R
\name{offset_sacling}
\alias{offset_sacling}
\title{Suitably generates offset matrix for the multivariate regression problem}
\usage{
offset_sacling(Y, ofset)
}
\arguments{
\item{Y}{outcome matrix}

\item{ofset}{offset matrix or microbiome data analysis specific scaling: common sum scaling = CSS (default), total sum scaling = TSS, median-ratio scaling = MRS, centered-log-ratio scaling  = CLR}
}
\description{
Suitably generates offset matrix for the multivariate regression problem
}
