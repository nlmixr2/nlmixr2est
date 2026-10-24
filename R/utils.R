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

#' A value C's `ISNA()` takes as missing: `NA`, but not `NaN`
#' @param x numeric vector
#' @return logical, one per element
#' @noRd
.nmIsNA <- function(x) {
  is.na(x) & !is.nan(x)
}

#' A data column as a double, as FOCEi reads DV and TIME (so "." is missing)
#' @param x data column
#' @return double vector
#' @noRd
.nmAsDouble <- function(x) {
  suppressWarnings(as.double(x))
}

#' Which records turn a compartment off: a negative or "-name" CMT
#' @param cmt the CMT column
#' @return logical, one per row
#' @noRd
.nmCmtOff <- function(cmt) {
  if (is.numeric(cmt)) {
    return(!is.na(cmt) & cmt < 0)
  }
  .c <- as.character(cmt)
  !is.na(.c) & startsWith(.c, "-")
}

#' Stop on event columns `etTrans()` rejects, with its message
#'
#' EVID, MDV and AMT must be stored as numbers: double, integer or logical, a
#' factor too (its codes), as `etTrans()` checks the storage type and not the
#' values.  No TIME may be infinite, and no ID missing.  A character TIME or DV
#' is read as a number, as FOCEi reads it.  The checks run in `etTrans()`'s order
#' and stop with its messages, after FOCEi's defaults: an ID it reads through
#' factor levels (any type but double) is checked first, and a double ID record
#' by record with TIME, after the event columns.  A logical ID, which `etTrans()`
#' rejects with a type-conversion error, stops first with a message naming it.
#' Without an EVID column but with AMT, the MDV column stands for the EVID and is
#' named "missing DV"; without an ID column the subject is 1.
#' @param d data frame with upper-case names
#' @return `d`, invisibly
#' @noRd
.nmAssertEventColumns <- function(d) {
  .id <- d[["ID"]]
  if (!is.null(.id) && !is.double(.id)) {
    .nmAssertIds(d)
  }
  .num <- c("double", "integer", "logical")
  if (!is.null(d[["EVID"]]) && !(typeof(d[["EVID"]]) %in% .num)) {
    stop("event id ('evid') needs to be an integer", call. = FALSE)
  }
  if (!is.null(d[["MDV"]]) && !(typeof(d[["MDV"]]) %in% .num)) {
    if (is.null(d[["EVID"]]) && !is.null(d[["AMT"]])) {
      stop("missing DV ('mdv') needs to be an integer", call. = FALSE)
    }
    stop("missing dependent variable ('mdv') needs to be an integer", call. = FALSE)
  }
  if (!is.null(d[["AMT"]]) && !(typeof(d[["AMT"]]) %in% .num)) {
    stop("amount ('amt') needs to be a number", call. = FALSE)
  }
  # record by record: the first infinite TIME or missing double ID (with no TIME
  # column, the first missing double ID)
  .time <- d[["TIME"]]
  .bad <- if (is.null(.time)) rep(FALSE, nrow(d)) else is.infinite(.nmAsDouble(.time))
  if (is.double(.id)) {
    .bad <- .bad | is.na(.id)
  }
  .row <- which(.bad)[1L]
  if (!is.na(.row)) {
    if (is.double(.id) && is.na(.id[.row])) {
      .nmAssertIds(d)
    }
    stop(
      "infinite times are not allowed (id: ",
      if (is.null(.id)) 1L else .id[.row],
      ", row: ",
      .row,
      ")",
      call. = FALSE
    )
  }
  invisible(d)
}

#' Which records of an event table are observations
#'
#' The records `etTrans()` keeps as observations (EVID 0), so data that has not
#' been through it agrees with what a solve returns, in every case
#' test-utils.R compares: a TIME; an EVID of 0 or missing; an MDV other than 1;
#' a CMT that does not turn a compartment off; and a DV -- or, with no DV, a
#' nonzero AMT, which `etTrans()` keeps as an observation (dropping the dose).
#' DV, TIME and AMT are read with `as.double()`, so a DV of "." is missing, and
#' EVID and MDV are truncated to integers, as `etTrans()` reads them; like its
#' `ISNA()`, a `NaN` is not missing.
#' @param d data frame with upper-case `TIME`, `EVID` and `DV`, and optionally
#'   `MDV`, `CMT` and `AMT`; [.nmWithEvid()] gives data without an EVID column one
#' @return logical, one per row
#' @noRd
.nmObsRows <- function(d) {
  if (is.null(d[["DV"]]) || is.null(d[["EVID"]]) || is.null(d[["TIME"]])) {
    stop("the data needs DV, EVID and TIME columns to find its observations", call. = FALSE)
  }
  .nmAssertEventColumns(d)
  .evid <- suppressWarnings(as.integer(.nmAsDouble(d[["EVID"]])))
  .obs <- !.nmIsNA(.nmAsDouble(d[["TIME"]])) & (is.na(.evid) | .evid == 0L)
  .mdv <- d[["MDV"]]
  if (!is.null(.mdv)) {
    .obs <- .obs & !(suppressWarnings(as.integer(.nmAsDouble(.mdv))) %in% 1L)
  }
  .cmt <- d[["CMT"]]
  if (!is.null(.cmt)) {
    .obs <- .obs & !.nmCmtOff(.cmt)
  }
  .amt <- d[["AMT"]]
  .dose <- if (is.null(.amt)) {
    FALSE
  } else {
    .a <- .nmAsDouble(.amt)
    !.nmIsNA(.a) & !(.a %in% 0)
  }
  .obs & (!.nmIsNA(.nmAsDouble(d[["DV"]])) | .dose)
}

#' Data with the EVID column FOCEi and `etTrans()` give data without one
#'
#' With no AMT column either, every record is EVID 0, as `.foceiPreProcessData()`
#' sets it.  Otherwise `etTrans()` takes an MDV column as the EVID (a missing MDV
#' is a dose, and MDV 1 with a zero AMT is EVID 2) and no longer reads it as an
#' MDV, so it is dropped; without one a nonzero AMT is a dose.
#' @param d data frame with upper-case names
#' @return `d`, with an integer EVID column when it had none
#' @noRd
.nmWithEvid <- function(d) {
  .nmAssertEventColumns(d)
  if (!is.null(d[["EVID"]])) {
    return(d)
  }
  .amt <- d[["AMT"]]
  if (is.null(.amt)) {
    d$EVID <- rep(0L, nrow(d))
    return(d)
  }
  .amt <- .nmAsDouble(.amt)
  .mdv <- d[["MDV"]]
  if (is.null(.mdv)) {
    d$EVID <- ifelse(.nmIsNA(.amt) | .amt %in% 0, 0L, 1L)
    return(d)
  }
  .evid <- suppressWarnings(as.integer(.nmAsDouble(.mdv)))
  .evid[.evid %in% 1L & .amt %in% 0] <- 2L
  .evid[is.na(.evid)] <- 1L
  d$EVID <- .evid
  d$MDV <- NULL
  d
}

#' Stop on a logical ID, or a missing one
#'
#' @param d data frame with an upper-case `ID` column
#' @return `d`, invisibly
#' @noRd
.nmAssertIds <- function(d) {
  .id <- d[["ID"]]
  # etTrans() cannot read a logical ID at all; its own error is a type-conversion
  # one ("Not compatible with requested type"), so this one names the column
  if (is.logical(.id)) {
    stop("'id' cannot be logical", call. = FALSE)
  }
  if (anyNA(.id)) {
    # etTrans()'s messages: it reads a double ID record by record, any other
    # through its factor levels
    stop(if (is.double(.id)) "'id' cannot be 'NA'" else "'id' cannot have NA values", call. = FALSE)
  }
  invisible(d)
}

#' Stop on data without an ID column, then where FOCEi stops, in its order
#'
#' The ID column, then [.nmAssertEventColumns()], which checks the IDs where
#' `etTrans()` does.  The VAE and VI preps need subjects, so data without an ID
#' column is an error there rather than FOCEi's single subject.
#' @param d data frame with upper-case names
#' @return `d`, invisibly
#' @noRd
.nmAssertIdData <- function(d) {
  if (is.null(d[["ID"]])) {
    stop("the data needs an ID column", call. = FALSE)
  }
  .nmAssertEventColumns(d)
}

#' Drop the subjects with no observations, as the FOCEi data step does
#'
#' Mirrors `.foceiPreProcessData()`: the subjects are named in a warning, and data
#' with no observation at all is kept whole.
#' @param d data frame with upper-case `ID`, `TIME`, `EVID` and `DV`
#' @return `d` without those subjects' records
#' @noRd
.nmDropNoObsSubjects <- function(d) {
  .nmAssertIdData(d)
  .obsId <- unique(d[["ID"]][.nmObsRows(d)])
  .drop <- setdiff(unique(d[["ID"]]), .obsId)
  if (length(.drop) == 0L || length(.obsId) == 0L) {
    return(d)
  }
  warning("IDs without observations dropped: ", paste(.drop, collapse = " "), call. = FALSE)
  d[d[["ID"]] %in% .obsId, , drop = FALSE]
}

#' Each subject's records and observations
#'
#' A subject's records are put in time order within each reset (EVID 3 or 4)
#' episode (ties keep their order), which is the order a solve returns its
#' observations in.
#' @param d data frame with upper-case `ID`, `TIME`, `EVID` and `DV`, and optionally
#'   `MDV`, `CENS` and `LIMIT`
#' @param ids the subjects wanted, in order, as values of `d$ID`
#' @return one element per `ids`: `NULL` for a subject with no records, else
#'   `list(ev, times, y, n, cens, limit)` -- its records, then the time, DV, count,
#'   censoring code (0 without `CENS`) and limit (`NA` without `LIMIT`) of its
#'   observations ([.nmObsRows()])
#' @noRd
.nmSubjectObs <- function(d, ids) {
  .nmAssertIdData(d)
  if (is.null(d[["TIME"]])) {
    stop("the data needs a TIME column", call. = FALSE)
  }
  .byId <- split(d, as.character(d[["ID"]]))
  lapply(as.character(ids), function(.id) {
    .di <- .byId[[.id]]
    if (is.null(.di) || nrow(.di) == 0L) {
      return(NULL)
    }
    # a reset (EVID 3 or 4) starts an episode that may reuse earlier clock times
    .episode <- cumsum(suppressWarnings(as.integer(.nmAsDouble(.di[["EVID"]]))) %in% c(3L, 4L))
    .di <- .di[order(.episode, suppressWarnings(as.double(.di[["TIME"]]))), , drop = FALSE]
    .obs <- .di[.nmObsRows(.di), , drop = FALSE]
    .n <- nrow(.obs)
    list(
      ev = .di,
      times = as.double(.obs[["TIME"]]),
      y = suppressWarnings(as.double(.obs[["DV"]])),
      n = .n,
      cens = if (is.null(.obs[["CENS"]])) integer(.n) else as.integer(.obs[["CENS"]]),
      limit = if (is.null(.obs[["LIMIT"]])) rep(NA_real_, .n) else as.numeric(.obs[["LIMIT"]])
    )
  })
}

#' Message suffix for a method that requires random effects
#'
#' Used as the `extra` argument to `rxode2::assertRxUiMixedOnly()` so a
#' single-subject / fixed-effect ("N of 1") model gets an actionable error
#' pointing to the methods that can fit models with no random effects.
#'
#' @param est estimation method name
#' @return character message suffix
#' @author Matthew L. Fidler
#' @noRd
.noRandomEffectMsg <- function(est) {
  paste0(
    " for the estimation routine '",
    est,
    "'; a model with no random effects (for example single-subject or 'N of 1' data) can be ",
    "fit with 'focei', 'foce', or a population method such as 'nlminb', 'bobyqa' or 'nls'"
  )
}

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
    .ret <- .lst[[arg]]
  }
  .ret
}


# ####################################################################### #
#
## Utilities for building nlmixr2
#
# ####################################################################### #

.rbindParHistory <- function(p1, p2) {
  .ret <- try(rbind(p1, p2), silent = TRUE)
  if (inherits(.ret, "try-error")) {
    warning("parameter history may be incomplete")
    .ret <- p2
  }
  .ret
}

refresh <- function() {
  ## nocov start
  source(devtools::package_file("build/refresh.R"))
  ## nocov end
}

nsis <- function() {
  ## build installer...
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
#' @param collectErr When \code{TRUE}, errors raised during evaluation
#'     of \code{expr} are recorded in addition to warnings.  A calling
#'     handler captures every error message as it is signalled, then
#'     lets the condition continue to propagate so that inner
#'     \code{try()}/\code{tryCatch()} blocks in \code{expr} keep
#'     working as usual.  An outer \code{tryCatch()} catches errors
#'     that escape all inner handling, so the call always returns
#'     instead of stopping.  If \code{expr} evaluated to a result
#'     (i.e. no error escaped), the recorded errors were caught by
#'     inner handlers and are discarded.  If an error did escape,
#'     every message observed along the error chain (including
#'     follow-up errors raised by \code{on.exit} handlers, see issue
#'     607) is returned in the \code{error} element of the result list
#'     (when \code{lst = TRUE}) or re-raised joined by newlines (when
#'     \code{lst = FALSE}).  When \code{FALSE} (the default) errors
#'     propagate normally and \code{es} is always \code{NULL}.  This
#'     is used by \code{nlmixr2Est0()} so that all errors from a
#'     failed estimation run are reported together rather than only
#'     the last one.
#'
#' @return The value of the expression, or when \code{lst = TRUE} a
#'     list of the form \code{list(object, warning = ws, error = es)}
#'     where \code{ws} and \code{es} are character vectors of unique
#'     warning and error messages (\code{es} is always \code{NULL}
#'     when \code{collectErr = FALSE}, and is also \code{NULL} when
#'     the expression evaluated successfully under \code{collectErr =
#'     TRUE}).
#'
#' @author Matthew L. Fidler
#' @keywords internal
#' @export
.collectWarn <- function(expr, lst = FALSE, collectErr = FALSE) {
  ws <- NULL
  es <- NULL
  if (getOption("nlmixr2.collectWarnings", TRUE)) {
    this.env <- environment()
    if (collectErr) {
      errored <- FALSE
      ret <- tryCatch(
        suppressWarnings(withCallingHandlers(
          expr,
          warning = function(w) {
            assign("ws", unique(c(w$message, ws)), this.env)
          },
          error = function(e) {
            assign("es", unique(c(e$message, es)), this.env)
          }
        )),
        error = function(e) {
          assign("errored", TRUE, this.env)
          NULL
        }
      )
      if (!errored) {
        es <- NULL
      }
    } else {
      ret <-
        suppressWarnings(withCallingHandlers(
          expr,
          warning = function(w) {
            assign("ws", unique(c(w$message, ws)), this.env)
          }
        ))
    }
  } else {
    ret <- force(expr)
  }
  if (lst) {
    ret <- list(ret, warning = ws, error = es)
  } else {
    for (w in ws) {
      warning(w)
    }
    if (length(es) > 0) {
      stop(paste(es, collapse = "\n"))
    }
  }
  ret
}
# #########################################################################

# nlmixr2Print() -----------------------------------------------------------
#' Print x using the message facility
#'
#' Captures print() output and routes it through message() so
#' suppressMessages() works on print functions.
#' @param x object to print
#' @return Nothing, called for its side effects
#' @param ... Other things output
#' @author Matthew L. Fidler
#' @export
#' @keywords internal
nlmixr2Print <- function(x, ...) {
  this.env <- environment()
  message(
    invisible(paste(
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
      warning("improper argument for 'rho'", call. = FALSE)
    }
    rho <- environment(fr)
  }
  step <- -.2 * start

  con <- list(maxeval = 999, reltol = 1e-6, rcoeff = 1., ecoeff = 2., ccoeff = .5, trace = FALSE) # nolint
  nmsC <- names(con)
  con[(namc <- names(control))] <- control
  if (length(noNms <- namc[!namc %in% nmsC])) {
    warning("unknown names in control: ", paste(noNms, collapse = ", "))
  }

  .Call(
    neldermead_wrap,
    fr,
    rho,
    length(start),
    start,
    step,
    as.integer(con$maxeval),
    con$reltol,
    con$rcoeff,
    con$ecoeff,
    con$ccoeff,
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
#' Makes C-level `REprintf` respect `suppressMessages()`, like R messages.
#'
#' @return Nothing
#' @keywords internal
#' @author Matthew Fidler
#' @export
#' @examples
#'
#' # Called automatically by nlmixr2(); other packages can call it too.
nmSuppressMsg <- function() {
  if (requireNamespace("knitr", quietly = TRUE)) {
    if (!is.null(knitr::opts_knit$get("rmarkdown.pandoc.to"))) {
      return(invisible(NULL))
    } else {
      .Call(
        `_nlmixr2est_setSilentErr`,
        as.integer(length(capture.output(message(" "), type = "message")) == 0L),
        PACKAGE = "nlmixr2est"
      )
    }
  } else {
    .Call(
      `_nlmixr2est_setSilentErr`,
      as.integer(length(capture.output(message(" "), type = "message")) == 0L),
      PACKAGE = "nlmixr2est"
    )
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
#' @param ensureSymmetry logical; symmetrizes `x` via
#'   \code{\link[Matrix]{symmpart}} unless already symmetric.
#'   \emph{Beware}: setting \code{FALSE} for asymmetric input is
#'   typically nonsense.
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
nmNearPD <- function(
  x,
  keepDiag = FALSE,
  do2eigen = TRUE,
  doDykstra = TRUE,
  only.values = FALSE,
  ensureSymmetry = !isSymmetric(x),
  eig.tol = 1e-6,
  conv.tol = 1e-7,
  posd.tol = 1e-8,
  maxit = 100L,
  trace = FALSE # nolint
) {
  if (ensureSymmetry) {
    x <- 0.5 * (t(x) + x)
  }
  .Call(
    `_nlmixr2est_nmNearPD_`,
    x,
    keepDiag,
    do2eigen,
    doDykstra,
    only.values,
    eig.tol,
    conv.tol,
    posd.tol,
    maxit,
    trace # nolint
  )
}

.sampleOmega <- function(omega) {
  rxode2::rxRmvn(1, sigma = omega)
}

# The fit's rxControl(cores=) for rxOptExpr(parallel=): chunked expression
# optimization then parallelizes with the same thread setting the solves use
# (0 keeps rxControl(cores=)'s meaning, the rxode2 thread setting).
.optExprCores <- function(ui) {
  .cores <- tryCatch(as.integer(rxode2::rxGetControl(ui, "rxControl", rxode2::rxControl())$cores), error = function(e) {
    0L
  })
  if (!length(.cores) || is.na(.cores)) 0L else .cores
}
