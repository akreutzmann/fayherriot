% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/estim_sigmau2.R
\name{wrapper_estsigmau2}
\alias{wrapper_estsigmau2}
\title{Wrapper function for the estmation of sigmau2}
\usage{
wrapper_estsigmau2(framework, method, interval)
}
\arguments{
\item{interval}{interval for the algorithm.}

\item{vardir}{direct variance.}

\item{precision}{precision criteria for the estimation of sigmau2.}

\item{maxiter}{maximum of iterations for the estimation of sigmau2.}

\item{direct}{direct estimator.}

\item{x}{matrix with explanatory variables.}

\item{areanumber}{number of domains.}
}
\value{
estimated sigmau2.
}
\description{
This function wraps the different estimation methods for sigmau2.
}
\keyword{internal}
