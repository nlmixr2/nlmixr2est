% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/resid.R
\name{addTable}
\alias{addTable}
\title{Add table information to nlmixr2 fit object without tables}
\usage{
addTable(
  object,
  updateObject = FALSE,
  data = object$dataSav,
  thetaEtaParameters = object$foceiThetaEtaParameters,
  table = tableControl(),
  keep = NULL,
  drop = NULL,
  envir = parent.frame(1)
)
}
\arguments{
\item{object}{nlmixr2 family of objects}

\item{updateObject}{Update the object (default FALSE)}

\item{data}{Saved data from}

\item{thetaEtaParameters}{Internal theta/eta parameters}

\item{table}{a `tableControl()` list of options}

\item{keep}{Character Vector of items to keep}

\item{drop}{Character Vector of items to drop or NULL}

\item{envir}{Environment to search for updating}
}
\value{
Fit with table information attached
}
\description{
Add table information to nlmixr2 fit object without tables
}
\examples{

\donttest{

one.cmt <- function() {
  ini({
    ## You may label each parameter with a comment
    tka <- 0.45 # Log Ka
    tcl <- log(c(0, 2.7, 100)) # Log Cl
    ## This works with interactive models
    ## You may also label the preceding line with label("label text")
    tv <- 3.45; label("log V")
    ## the label("Label name") works with all models
    eta.ka ~ 0.6
    eta.cl ~ 0.3
    eta.v ~ 0.1
    add.sd <- 0.7
  })
  model({
    ka <- exp(tka + eta.ka)
    cl <- exp(tcl + eta.cl)
    v <- exp(tv + eta.v)
    linCmt() ~ add(add.sd)
  })
}

# run without tables step
f <- nlmixr2(one.cmt, theo_sd, "saem", control=list(calcTables=FALSE))

print(f)

# Now add the tables

f <- addTable(f)

print(f)

}
}
\author{
Matthew Fidler
}
