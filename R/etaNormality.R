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
  .ret$meanP <- stats::t.test(x, mu = 0)$p.value
  # shapiro.test() only accepts 3 to 5000 values
  if (.n <= 5000L) {
    .sw <- stats::shapiro.test(x)
    .ret$shapiroW <- unname(.sw$statistic)
    .ret$shapiroP <- .sw$p.value
  }
  .ret
}

#' Normality tests for each ETA column of a data frame
#'
#' @param eta data frame of ETAs (an `ID`, `mixnum` or occasion column is ignored)
#' @param skip names of ETAs whose statistics are `NA` (fixed or zero omega)
#' @param drop names of non-ETA columns to ignore
#' @return data frame, one row per ETA column (see `.etaNormalityRow()`)
#' @author Matthew L. Fidler
#' @noRd
.etaNormalityDf <- function(eta, skip = character(0), drop = character(0)) {
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
      .etaNormalityRow(.c, eta[[.c]], skip = .c %in% skip)
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

#' @rdname nmObjGet
#' @export
nmObjGet.etaNormality <- function(x, ...) {
  .env <- x[[1]]
  .ranef <- NULL
  if (exists("ranef", envir = .env, inherits = FALSE)) {
    .ranef <- get("ranef", envir = .env, inherits = FALSE)
  }
  .ui <- NULL
  if (exists("ui", envir = .env, inherits = FALSE)) {
    .ui <- rxode2::rxUiDecompress(get("ui", envir = .env, inherits = FALSE))
  }
  .skip <- .etaNormalityFixed(.ui$iniDf)
  .ret <- .etaNormalityDf(.ranef, skip = .skip)
  # occasion etas, pooled over every subject and occasion
  if (exists("iov", envir = .env, inherits = FALSE)) {
    .iov <- get("iov", envir = .env, inherits = FALSE)
    if (is.list(.iov) && !is.data.frame(.iov)) {
      for (.occ in names(.iov)) {
        .ret <- rbind(.ret, .etaNormalityDf(.iov[[.occ]], skip = .skip, drop = .occ))
      }
    }
  }
  if (!is.null(.ret)) {
    rownames(.ret) <- NULL
  }
  .ret
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
