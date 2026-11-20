#' One ETA's mean, sd and zero-mean t-test p-value from `$shrink`
#'
#' @param shrink the fit's `shrink`: a data frame, or (IOV) a list of them
#' @param eta name of the ETA
#' @return named numeric vector `mean`, `sd`, `meanP` (`NA` when absent)
#' @author Matthew L. Fidler
#' @noRd
.etaNormalityShrink <- function(shrink, eta) {
  .ret <- c(mean = NA_real_, sd = NA_real_, meanP = NA_real_)
  .tabs <- if (is.data.frame(shrink)) list(shrink) else shrink
  for (.t in .tabs) {
    if (is.data.frame(.t) && eta %in% names(.t)) {
      .v <- .t[[eta]]
      names(.v) <- rownames(.t)
      .ret[] <- unname(.v[c("mean", "sd", "p-value")])
      break
    }
  }
  .ret
}

#' Shapiro-Wilk test of one ETA's estimates
#'
#' @param x numeric vector of the ETA's empirical Bayes estimates
#' @return named numeric vector `shapiroW`, `shapiroP`; `NA` outside the 3 to
#'   5000 values `shapiro.test()` accepts, or for constant estimates
#' @author Matthew L. Fidler
#' @noRd
.etaNormalityShapiro <- function(x) {
  .ret <- c(shapiroW = NA_real_, shapiroP = NA_real_)
  if (length(x) < 3L || length(x) > 5000L) {
    return(.ret)
  }
  .sw <- tryCatch(stats::shapiro.test(x), error = function(e) NULL)
  if (!is.null(.sw)) {
    .ret[] <- c(unname(.sw$statistic), .sw$p.value)
  }
  .ret
}

#' Normality row for one ETA
#'
#' @param eta name of the ETA
#' @param x numeric vector of the ETA's empirical Bayes estimates
#' @param shrink the fit's `shrink` (source of `mean`, `sd` and `meanP`)
#' @param skip when `TRUE` (a fixed or zero-variance ETA), every statistic is `NA`
#' @return one-row data frame with columns `eta`, `n`, `mean`, `sd`, `meanP`,
#'   `shapiroW` and `shapiroP`
#' @author Matthew L. Fidler
#' @noRd
.etaNormalityRow <- function(eta, x, shrink = NULL, skip = FALSE) {
  x <- x[is.finite(x)]
  .s <- c(mean = NA_real_, sd = NA_real_, meanP = NA_real_, shapiroW = NA_real_, shapiroP = NA_real_)
  if (!skip) {
    .s <- c(.etaNormalityShrink(shrink, eta), .etaNormalityShapiro(x))
  }
  data.frame(eta = eta, n = length(x), as.list(.s), stringsAsFactors = FALSE)
}

#' Normality rows for each ETA column of a data frame
#'
#' @param eta data frame of ETAs (an `ID`, `mixnum` or occasion column is ignored)
#' @param shrink the fit's `shrink`
#' @param skip names of ETAs whose statistics are `NA` (fixed or zero omega)
#' @param drop names of non-ETA columns to ignore
#' @param rows named list of the rows to use for an ETA (default every row)
#' @return data frame, one row per ETA column (see `.etaNormalityRow()`)
#' @author Matthew L. Fidler
#' @noRd
.etaNormalityDf <- function(eta, shrink = NULL, skip = character(0), drop = character(0), rows = list()) {
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
      .etaNormalityRow(.c, .x, shrink = shrink, skip = .c %in% skip)
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
#' @param fit fit environment, or a list with `ui`, `ranef`, `shrink` and
#'   optionally `iov` and `mixNum`
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
  .shrink <- .get("shrink")
  .rows <- .etaNormalityMixRows(.ui, .ranef, .mixNum)
  .ret <- .etaNormalityDf(.ranef, shrink = .shrink, skip = .skip, rows = .rows)
  # occasion etas: Shapiro-Wilk pooled over every subject and occasion
  .iov <- .get("iov")
  if (is.list(.iov) && !is.data.frame(.iov)) {
    for (.occ in names(.iov)) {
      .d <- .iov[[.occ]]
      .ret <- rbind(
        .ret,
        .etaNormalityDf(
          .d,
          shrink = .shrink,
          skip = .skip,
          drop = .occ,
          rows = .etaNormalityMixRows(.ui, .d, .mixNum)
        )
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
attr(nmObjGet.etaNormality, "desc") <- "ETA normality (Shapiro-Wilk) test with the $shrink mean test"
attr(nmObjGet.etaNormality, "rstudio") <- data.frame(
  eta = "eta.cl",
  n = 3L,
  mean = 0,
  sd = 1,
  meanP = 1,
  shapiroW = 1,
  shapiroP = 1
)
