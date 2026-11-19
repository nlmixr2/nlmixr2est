#' Normality and zero-mean tests for one ETA's estimates
#'
#' @param eta name of the ETA
#' @param x numeric vector of the ETA's empirical Bayes estimates
#' @param skip when `TRUE` (a fixed or zero-variance ETA), every statistic is `NA`
#' @return one-row data frame with columns `eta`, `n`, `mean`, `sd`, `meanP`,
#'   `shapiroW` and `shapiroP`
#' @author Matthew L. Fidler
#' @noRd
.etaNormalityRow <- function(eta, x, skip = FALSE) {
  x <- x[is.finite(x)]
  .n <- length(x)
  .ret <- data.frame(
    eta = eta,
    n = .n,
    mean = NA_real_,
    sd = NA_real_,
    meanP = NA_real_,
    shapiroW = NA_real_,
    shapiroP = NA_real_,
    stringsAsFactors = FALSE
  )
  if (skip || .n < 3L) {
    return(.ret)
  }
  .ret$mean <- mean(x)
  .ret$sd <- stats::sd(x)
  # constant estimates (all shrunk to one value) have no test
  if (!is.finite(.ret$sd) || .ret$sd <= 0) {
    return(.ret)
  }
  # t.test() errors on essentially constant data
  .ret$meanP <- tryCatch(stats::t.test(x, mu = 0)$p.value, error = function(e) NA_real_)
  # shapiro.test() only accepts 3 to 5000 values
  if (.n <= 5000L) {
    .sw <- tryCatch(stats::shapiro.test(x), error = function(e) NULL)
    if (!is.null(.sw)) {
      .ret$shapiroW <- unname(.sw$statistic)
      .ret$shapiroP <- .sw$p.value
    }
  }
  .ret
}

#' Normality tests for each ETA column of a data frame
#'
#' @param eta data frame of ETAs (an `ID`, `mixnum` or occasion column is ignored)
#' @param skip names of ETAs whose statistics are `NA` (fixed or zero omega)
#' @param drop names of non-ETA columns to ignore
#' @param rows named list of the rows to use for an ETA (default every row)
#' @return data frame, one row per ETA column (see `.etaNormalityRow()`)
#' @author Matthew L. Fidler
#' @noRd
.etaNormalityDf <- function(eta, skip = character(0), drop = character(0), rows = list()) {
  if (!is.data.frame(eta)) {
    return(NULL)
  }
  eta <- .nmDropNonEtaCols(eta)
  .cols <- names(eta)[vapply(eta, is.numeric, logical(1))]
  .cols <- setdiff(.cols, drop)
  if (length(.cols) == 0L) {
    return(NULL)
  }
  do.call(
    rbind,
    lapply(.cols, function(.c) {
      .x <- eta[[.c]]
      if (!is.null(rows[[.c]])) {
        .x <- .x[rows[[.c]]]
      }
      .etaNormalityRow(.c, .x, skip = .c %in% skip)
    })
  )
}

#' Names of the ETAs that are fixed or have a zero variance
#'
#' @param iniDf the fit's final `iniDf`
#' @return character vector of ETA names
#' @author Matthew L. Fidler
#' @noRd
.etaNormalityFixed <- function(iniDf) {
  if (!is.data.frame(iniDf)) {
    return(character(0))
  }
  .w <- which(!is.na(iniDf$neta1) & iniDf$neta1 == iniDf$neta2)
  .w <- .w[iniDf$fix[.w] | iniDf$est[.w] == 0]
  iniDf$name[.w]
}

#' Subjects whose estimates a component-owned mixture ETA uses
#'
#' An ETA used only in one mixture component is 0 for the other components'
#' subjects, so only its own component's subjects are tested.
#'
#' @param ui the fit's final ui
#' @param ranef data frame of ETAs with an `ID` column
#' @param mixNum the subjects' components (`ID` and `mixnum` columns)
#' @return named list of logical row indices, one per owned ETA
#' @author Matthew L. Fidler
#' @noRd
.etaNormalityMixRows <- function(ui, ranef, mixNum) {
  if (is.null(ui) || !is.data.frame(ranef) || !is.data.frame(mixNum)) {
    return(list())
  }
  .owner <- .mixEtaOwner(ui)
  .owner <- .owner[.owner > 0L & names(.owner) %in% names(ranef)]
  .mn <- mixNum$mixnum[match(ranef$ID, mixNum$ID)]
  lapply(.owner, function(.k) !is.na(.mn) & .mn == .k)
}

#' ETA normality table for a fit
#'
#' @param fit fit environment, or a list with `ui`, `ranef` and optionally
#'   `iov` and `mixNum`
#' @return data frame, one row per ETA (see `.etaNormalityRow()`)
#' @author Matthew L. Fidler
#' @noRd
.etaNormalityFit <- function(fit) {
  .get <- function(n) {
    if (is.environment(fit)) {
      return(get0(n, envir = fit, inherits = FALSE))
    }
    fit[[n]]
  }
  .ranef <- .get("ranef")
  .ui <- .get("ui")
  if (!is.null(.ui)) {
    .ui <- rxode2::rxUiDecompress(.ui)
  }
  .skip <- .etaNormalityFixed(.ui$iniDf)
  .mixNum <- .get("mixNum")
  if (!is.data.frame(.mixNum) && is.data.frame(.ranef) && "mixnum" %in% names(.ranef)) {
    .mixNum <- .ranef[, c("ID", "mixnum"), drop = FALSE]
  }
  .rows <- .etaNormalityMixRows(.ui, .ranef, .mixNum)
  .ret <- .etaNormalityDf(.ranef, skip = .skip, rows = .rows)
  # occasion etas, pooled over every subject and occasion
  .iov <- .get("iov")
  if (is.list(.iov) && !is.data.frame(.iov)) {
    for (.occ in names(.iov)) {
      .d <- .iov[[.occ]]
      .ret <- rbind(
        .ret,
        .etaNormalityDf(.d, skip = .skip, drop = .occ, rows = .etaNormalityMixRows(.ui, .d, .mixNum))
      )
    }
  }
  # zero-variance etas are dropped from the fit's etas, so give them NA rows
  .miss <- setdiff(.skip, .ret$eta)
  if (length(.miss) > 0L) {
    .add <- do.call(rbind, lapply(.miss, .etaNormalityRow, x = numeric(0), skip = TRUE))
    .add$n <- NA_integer_
    .ret <- rbind(.ret, .add)
  }
  if (!is.null(.ret)) {
    rownames(.ret) <- NULL
  }
  .ret
}

#' @rdname nmObjGet
#' @export
nmObjGet.etaNormality <- function(x, ...) {
  .etaNormalityFit(x[[1]])
}
attr(nmObjGet.etaNormality, "desc") <- "ETA mean (t-test) and normality (Shapiro-Wilk) tests"
attr(nmObjGet.etaNormality, "rstudio") <- data.frame(
  eta = "eta.cl",
  n = 3L,
  mean = 0,
  sd = 1,
  meanP = 1,
  shapiroW = 1,
  shapiroP = 1
)
