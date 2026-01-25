## utils.R: population PK/PD modeling library
##
## Copyright (C) 2014 - 2016  Wenping Wang,
##  Portions (c) Matt Fidler and rest of the nlmixr2 team.
##  The authorship is always up to date on git
##
## This file is part of nlmixr2.
##
## nlmixr2 is free software: you can redistribute it and/or modify it
## under the terms of the GNU General Public License as published by
## the Free Software Foundation, either version 2 of the License, or
## (at your option) any later version.
##
## nlmixr2 is distributed in the hope that it will be useful, but
## WITHOUT ANY WARRANTY; without even the implied warranty of
## MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
## GNU General Public License for more details.
##
## You should have received a copy of the GNU General Public License
## along with nlmixr2.  If not, see <http:##www.gnu.org/licenses/>.

# Utilities for nlmixr2 ####################################################

#' Cox Box, Yeo Johnson and inverse transformation
#'
#' @param x data to transform
#' @param lambda Cox-box lambda parameter
#' @return Cox-Box Transformed Data
#' @author Matthew L. Fidler
#' @examples
#'
#' boxCox(1:3,1) ## Normal
#' iBoxCox(boxCox(1:3,1))
#'
#' boxCox(1:3,0) ## Log-Normal
#' iBoxCox(boxCox(1:3,0),0)
#'
#' boxCox(1:3,0.5) ## lambda=0.5
#' iBoxCox(boxCox(1:3,0.5),0.5)
#'
#' yeoJohnson(seq(-3,3),1) ## Normal
#' iYeoJohnson(yeoJohnson(seq(-3,3),1))
#'
#' yeoJohnson(seq(-3,3),0)
#' iYeoJohnson(yeoJohnson(seq(-3,3),0),0)
#' @export
boxCox <- function(x, lambda = 1) {
  checkmate::assert_numeric(x)
  checkmate::assert_numeric(lambda)
  .Call(`_nlmixr2est_boxCox_`, x, lambda, 0L)
}

#' @rdname boxCox
#' @export
iBoxCox <- function(x, lambda = 1) {
  .Call(`_nlmixr2est_iBoxCox_`, x, lambda, 0L)
}

#' @rdname boxCox
#' @export
yeoJohnson <- function(x, lambda = 1) {
  .Call(`_nlmixr2est_boxCox_`, x, lambda, 1L)
}

#' @rdname boxCox
#' @export
iYeoJohnson <- function(x, lambda = 1) {
  .Call(`_nlmixr2est_iBoxCox_`, x, lambda, 1L)
}


#' @importFrom utils capture.output
#' @noRd
.captureOutput <- function(expr, envir = parent.frame()) {
  eval(
    {
      .file <- rawConnection(raw(0L), open = "w")
      on.exit({
        if (!is.null(.file)) close(.file)
      })
      capture.output(expr, file = .file)
      .ret <- rawConnectionValue(.file)
      close(.file)
      .file <- NULL
      .ret <- rawToChar(.ret)
      return(.ret)
    },
    envir = envir,
    enclos = envir
  )
}

#' @export
`$.nlmixr2Gill83` <- function(obj, arg, exact = FALSE) {
  .ret <- obj[[arg]]
  if (is.null(.ret)) {
    .cls <- class(obj)
    .lst <- attr(.cls, ".nlmixr2Gill")
    return(.lst[[arg]])
  }
  .ret
}


# ####################################################################### #
#
## Utilities for building nlmixr2
#
# ####################################################################### #

.rbindParHistory <- function(p1, p2) {
  .ret <- try(rbind(p1, p2), silent=TRUE)
  if (inherits(.ret, "try-error")) {
    warning("parameter history may be incomplete")
    return(p2)
  }
  .ret
}

refresh <- function() {
  ## nocov start
  source(devtools::package_file("build/refresh.R"))
  ## nocov end
}

nsis <- function() { ## build installer...
  ## nocov start
  source(devtools::package_file("build/nsis.R"))
  ## nocov end
}
# ########################################################################

# .collectWarn --------------------------------------------------------
#' Collect warnings and just warn once.
#'
#' @param expr R expression
#'
#' @param lst When \code{TRUE} return a list with
#'     list(object,warnings) instead of issuing the warnings.
#'     Otherwise, when \code{FALSE} issue the warnings and return the
#'     object.
#'
#' @return The value of the expression or a list with the value of
#'     the expression and a list of warning messages
#'
#' @author Matthew L. Fidler
#'
#' @export
#'
#' @keywords internal
.collectWarn <- function(expr, lst = FALSE) {
  if (getOption("nlmixr2.collectWarnings", TRUE)) {
    ws <- c()
    this.env <- environment()
    ret <-
      suppressWarnings(withCallingHandlers(
        expr,
        warning = function(w) {
          assign("ws", unique(c(w$message, ws)), this.env)
        }
      ))
    if (lst) {
      list(ret, ws)
    } else {
      for (w in ws) {
        warning(w)
      }
      ret
    }
  } else {
    ret <- force(expr)
    if (lst) {
      return(list(ret, NULL))
    }
    ret
  }
}
# #########################################################################

# nlmixr2Print() -----------------------------------------------------------
#' Print x using the message facility
#'
#' This allows the suppressMessages to work on print functions.  This
#' captures the output function sends it through the message routine.
#'
#' catpureOutput was used since it is much faster than the internal
#' capture.output see https://www.r-bloggers.com/performance-captureoutput-is-much-faster-than-capture-output/
#' @param x object to print
#' @return Nothing, called for its side effects
#' @param ... Other things output
#' @author Matthew L. Fidler
#' @export
#' @keywords internal
nlmixr2Print <- function(x, ...) {
  this.env <- environment()
  message(invisible(paste(
    .captureOutput(assign("x", print(x, ...), this.env)),
    collapse = "\n"
  )),
  appendLF = TRUE
  )
  invisible(x)
}
# #########################################################################

# cholSE() ----------------------------------------------------------------
#' Generalized Cholesky Matrix Decomposition
#'
#'  Performs a (modified) Cholesky factorization of the form
#'
#'   t(P) \%*\% A \%*\% P  + E = t(R) \%*\% R
#'
#'  As detailed in Schnabel/Eskow (1990)
#'
#' @param matrix Matrix to be Factorized.
#' @param tol Tolerance; Algorithm suggests (.Machine$double.eps) ^ (1 / 3), default
#' @return Generalized Cholesky decomposed matrix.
#' @author Matthew L. Fidler (translation), Johannes Pfeifer, Robert
#'     B. Schnabel and Elizabeth Eskow
#'
#' @references
#'
#' matlab source: http://www.dynare.org/dynare-matlab-m2html/matlab/chol_SE.html; Slightly different return values
#'
#' Robert B. Schnabel and Elizabeth
#' Eskow. 1990. "A New Modified Cholesky Factorization," SIAM Journal
#' of Scientific Statistical Computing, 11, 6: 1136-58.
#'
#' Elizabeth Eskow and Robert B. Schnabel
#' 1991. "Algorithm 695 - Software for a New Modified Cholesky Factorization,"
#' ACM Transactions on Mathematical Software, Vol 17, No 3: 306-312
#'
#' @note
#'
#' This version does not pivot or return the E matrix
#'
#' @export
cholSE <- function(matrix, tol = (.Machine$double.eps)^(1 / 3)) {
  .Call(`_nlmixr2est_cholSE_`, matrix, tol)
}
# #########################################################################

.setRoot <- function() {
  setwd("c:/")
}

#' Generate a data.frame using the R4.0 convention
#'
#' @param ... Passed to \code{base::data.frame()} or
#'   \code{base::as.data.frame()}
#' @param stringsAsFactors Captured so that it can be ignored and always set to
#'   \code{FALSE}
#' @return A data.frame with strings not converted to factors
#' @noRd
.data.frame <- function(..., stringsAsFactors = FALSE) {
  base::data.frame(..., stringsAsFactors = FALSE)
}
.as.data.frame <- function(..., stringsAsFactors = FALSE) {
  base::as.data.frame(..., stringsAsFactors = FALSE)
}

#' Nelder-Mead simplex search
#'
#' @param start initials
#' @param fr objective function
#' @param rho evaluation environment
#' @param control additional optimization options
#' @return a list of ...
#' @export
nmsimplex <- function(start, fr, rho = NULL, control = list()) {
  if (!is.environment(rho)) {
    if (!is.null(rho)) {
      warning("improper argument for 'rho'", call.=FALSE)
    }
    rho <- environment(fr)
  }
  step <- -.2 * start

  con <- list(maxeval = 999, reltol = 1e-6, rcoeff = 1., ecoeff = 2., ccoeff = .5, trace = FALSE) # nolint
  nmsC <- names(con)
  con[(namc <- names(control))] <- control
  if (length(noNms <- namc[!namc %fin% nmsC])) {
    warning("unknown names in control: ", paste(noNms, collapse = ", "))
  }

  .Call(neldermead_wrap, fr, rho, length(start), start, step,
    as.integer(con$maxeval), con$reltol, con$rcoeff, con$ecoeff, con$ccoeff,
    as.integer(con$trace), # nolint
    PACKAGE = "nlmixr2est"
  )
}

#' Change a character expression into a quoted name
#'
#' @param chr Character expression
#' @return Quote name
#' @author Matthew L. Fidler
#' @noRd
.enQuote <- function(chr) {
  eval(parse(text = paste0("quote(", chr, ")")))
}

#' Respect suppress messages for nlmixr2 C functions
#'
#' This turns on the silent REprintf in C when `suppressMessages()` is
#' turned on. This makes the `REprintf` act like `messages` in R,
#' they can be suppressed with `suppressMessages()`
#'
#' @return Nothing
#' @keywords internal
#' @author Matthew Fidler
#' @export
#' @examples
#'
#' # nmSupressMsg() is called with nlmixr2()
#'
#' # In nlmixr2, we use REprintf so that interrupted threads do not crash R
#' # if there is a user interrupt. This isn't captured by R's messages, but
#' # This interface allows the `suppressMessages()` to suppress the C printing
#' # as well
#'
#' # If you  want to suppress messages from nlmixr2 in other packages, you can use
#' # this function
nmSuppressMsg <- function() {
  if (requireNamespace("knitr", quietly = TRUE)) {
    if (!is.null(knitr::opts_knit$get("rmarkdown.pandoc.to"))) {
      return(invisible(NULL))
    } else {
      .Call(`_nlmixr2est_setSilentErr`, as.integer(length(capture.output(message(" "), type = "message")) == 0L),
            PACKAGE="nlmixr2est")
    }
  } else {
    .Call(`_nlmixr2est_setSilentErr`, as.integer(length(capture.output(message(" "), type = "message")) == 0L),
          PACKAGE="nlmixr2est")
  }
  invisible(NULL)
}

#' @export
rxModelVarsS3.nlmixr2FitCore <- function(obj) {
  rxode2::rxModelVars(obj$ui)
}

#' @export
rxModelVarsS3.nlmixr2FitCoreSilent <- function(obj) {
  rxode2::rxModelVars(obj$ui)
}

#' C++ implementation of Matrix's nearPD
#'
#' With `ensureSymmetry` it makes sure it is symmetric by applying 0.5*(t(x) + x) before using nmNearPD
#'
#' @inherit Matrix::nearPD
#'
#' @param ensureSymmetry  logical; by default, \code{\link[Matrix]{symmpart}(x)}
#' is used whenever \code{isSymmetric(x)} is not true.  The user
#' can explicitly set this to \code{TRUE} or \code{FALSE}, saving the
#' symmetry test. \emph{Beware} however that setting it \code{FALSE}
#' for an \bold{a}symmetric input \code{x}, is typically nonsense!
#'
#' @return unlike the matrix package, this simply returns the nearest
#'   positive definite matrix
#'
#' @examples
#'
#' set.seed(27)
#' m <- matrix(round(rnorm(25),2), 5, 5)
#' m <- m + t(m)
#' diag(m) <- pmax(0, diag(m)) + 1
#' (m <- round(cov2cor(m), 2))
#'
#' near.m <- nmNearPD(m)
#' round(near.m, 2)
#' norm(m - near.m) # 1.102 / 1.08
#'
#' round(nmNearPD(m, only.values=TRUE), 9)
#'
#' ## A longer example, extended from Jens' original,
#' ## showing the effects of some of the options:
#'
#' pr <- matrix(c(1,     0.477, 0.644, 0.478, 0.651, 0.826,
#'                0.477, 1,     0.516, 0.233, 0.682, 0.75,
#'                0.644, 0.516, 1,     0.599, 0.581, 0.742,
#'                0.478, 0.233, 0.599, 1,     0.741, 0.8,
#'                0.651, 0.682, 0.581, 0.741, 1,     0.798,
#'                0.826, 0.75,  0.742, 0.8,   0.798, 1),
#'                nrow = 6, ncol = 6)
#'
#' nc  <- nmNearPD(pr)
#'
#' @export
nmNearPD <- function(x, keepDiag = FALSE, do2eigen = TRUE, doDykstra = TRUE, only.values = FALSE, ensureSymmetry=!isSymmetric(x), eig.tol = 1e-6, conv.tol = 1e-7, posd.tol = 1e-8, maxit = 100L,
                     trace = FALSE # nolint
                     ) {
  if (ensureSymmetry) {
    x <- 0.5 * (t(x) + x)
  }
  .Call(`_nlmixr2est_nmNearPD_`, x, keepDiag, do2eigen, doDykstra, only.values, eig.tol, conv.tol, posd.tol, maxit,
        trace # nolint
        )
}

.sampleOmega <- function(omega) {
  rxode2::rxRmvn(1, sigma=omega)
}
