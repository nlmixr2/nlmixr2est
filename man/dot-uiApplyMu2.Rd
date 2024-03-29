% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/mu2.R
\name{.uiApplyMu2}
\alias{.uiApplyMu2}
\title{This is an internal function for modifying the UI to apply mu2 referencing}
\usage{
.uiApplyMu2(env)
}
\arguments{
\item{env}{Environment needed for nlmixr2 fits}
}
\value{
Either the original model({}) block (if changed) or NULL if
  not changed
}
\description{
mu2 referencing is algebraic mu-referencing by converting to the
transformation to a single value in the original dataset, and
moving that around
}
\author{
Matthew L. Fidler
}
\keyword{internal}
