% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/plot-mfa.R
\name{plot_comp}
\alias{plot_comp}
\title{plot_comp}
\usage{
plot_comp(mfa, d = c(1, 2), ...)
}
\arguments{
\item{mfa}{an object of class \code{"mfa"}}

\item{d}{a vector of dimensions to display}

\item{\dots}{arguments to be passed to/from other methods}
}
\description{
Plots the compromise of an mfa object
}
\examples{
 \dontrun{
 # create an mfa object
 mfa1 <- mfa(data, sets=list(2:7,8:13,14:19)

 plot_comp(mfa1)
 }
}

