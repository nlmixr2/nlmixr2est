% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/focei.R
\name{nlmixr2CreateOutputFromUi}
\alias{nlmixr2CreateOutputFromUi}
\title{Create nlmixr output from the UI}
\usage{
nlmixr2CreateOutputFromUi(
  ui,
  data = NULL,
  control = NULL,
  table = NULL,
  env = NULL,
  est = "none"
)
}
\arguments{
\item{ui}{This is the UI that will be used for the translation}

\item{data}{This has the data}

\item{control}{focei control for data creation}

\item{table}{Table options}

\item{env}{Environment setup which needs the following:
- `$table` for table options
- `$origData` -- Original Data
- `$dataSav` -- Processed data from .foceiPreProcessData
- `$idLvl` -- Level information for ID factor added
- `$covLvl` -- Level information for items to convert to factor
- `$ui` for ui object
- `$fullTheta` Full theta information
- `$etaObf` data frame with ID, etas and OBJI
- `$cov` For covariance
- `$covMethod` for the method of calculating the covariance
- `$adjObf` Should the objective function value be adjusted
- `$objective` objective function value
- `$extra` Extra print information
- `$method` Estimation method (for printing)
- `$omega` Omega matrix
- `$theta` Is a theta data frame
- `$model` a list of model information for table generation.  Needs a `predOnly` model
- `$message` Message for display
- `$est` estimation method
- `$ofvType` (optional) tells the type of ofv is currently being use

There are some more details that need to be described here}

\item{est}{Estimation method}
}
\value{
nlmixr fit object
}
\description{
Create nlmixr output from the UI
}
\author{
Matthew L. Fidler
}
