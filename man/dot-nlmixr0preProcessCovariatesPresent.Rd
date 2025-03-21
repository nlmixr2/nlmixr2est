% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/preProcessCovariatesPresent.R
\name{.nlmixr0preProcessCovariatesPresent}
\alias{.nlmixr0preProcessCovariatesPresent}
\title{Preprocess Covariates needed (or other data items)}
\usage{
.nlmixr0preProcessCovariatesPresent(ui, est, data, control)
}
\arguments{
\item{ui}{rxode2 ui}

\item{est}{estimation method (all methods are shown by
`nlmixr2AllEst()`). Methods can be added for other tools}

\item{data}{nlmixr data}

\item{control}{The estimation control object.  These are expected
to be different for each type of estimation method}
}
\value{
list with the ui (possibly modified)
}
\description{
Preprocess Covariates needed (or other data items)
}
\author{
Matthew L. Fidler
}
