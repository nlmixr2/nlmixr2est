#' Shapiro-Wilk test of one ETA's estimates
#'
#' @param x numeric vector of the ETA's empirical Bayes estimates
#' @return named numeric vector `shapiroW`, `shapiroP`; `NA` outside the 3 to
#'   5000 values `shapiro.test()` accepts, or for constant estimates
#' @author Matthew L. Fidler
#' @noRd
.etaNormalityShapiro <- function(x) {
  .ret <- c(shapiroW = NA_real_, shapiroP = NA_real_)
  x <- x[is.finite(x)]
  if (length(x) < 3L || length(x) > 5000L) {
    return(.ret)
  }
  .sw <- tryCatch(stats::shapiro.test(x), error = function(e) NULL)
  if (!is.null(.sw)) {
    .ret[] <- c(unname(.sw$statistic), .sw$p.value)
  }
  .ret
}

# row names of the Shapiro-Wilk rows added to `$shrink`
.shrinkShapiroRows <- c("Shapiro-Wilk W", "Shapiro-Wilk p-value")

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

#' Each ETA column's estimates and the rows it uses
#'
#' @param eta data frame of ETAs (an `ID`, `mixnum` or `drop` column is ignored)
#' @param rows named list of the rows to use for an ETA (default every row)
#' @param drop names of non-ETA columns to ignore
#' @return named list of logical row indices, one per ETA column
#' @author Matthew L. Fidler
#' @noRd
.etaNormalityCols <- function(eta, rows = list(), drop = character(0)) {
  if (!is.data.frame(eta)) {
    return(list())
  }
  eta <- .nmDropNonEtaCols(eta)
  .cols <- setdiff(names(eta)[vapply(eta, is.numeric, logical(1))], drop)
  .ret <- lapply(.cols, function(.c) {
    if (is.null(rows[[.c]])) rep(TRUE, nrow(eta)) else rows[[.c]]
  })
  names(.ret) <- .cols
  .ret
}

#' The ETA estimates behind each `$shrink` column
#'
#' @param fit fit environment, or a list with `ranef` and optionally `ui`,
#'   `iov` and `mixNum`
#' @return list with `eta` (one vector per subject or occasion ETA, occasion
#'   ETAs pooled over occasions) and `occ` (one vector per `rx.<eta>.<occasion>`
#'   column of an IOV fit's `$shrink`)
#' @author Matthew L. Fidler
#' @noRd
.etaNormalityValues <- function(fit) {
  .get <- function(n) {
    if (is.environment(fit)) {
      return(get0(n, envir = fit, inherits = FALSE))
    }
    fit[[n]]
  }
  .ranef <- .get("ranef")
  .mixNum <- .get("mixNum")
  if (!is.data.frame(.mixNum) && is.data.frame(.ranef) && "mixnum" %in% names(.ranef)) {
    .mixNum <- .ranef[, c("ID", "mixnum"), drop = FALSE]
  }
  # only a mixture fit needs the ui (for its component-owned etas)
  .ui <- NULL
  if (is.data.frame(.mixNum)) {
    .ui <- .get("ui")
    if (!is.null(.ui)) {
      .ui <- rxode2::rxUiDecompress(.ui)
    }
  }
  .rows <- .etaNormalityCols(.ranef, rows = .etaNormalityMixRows(.ui, .ranef, .mixNum))
  .eta <- lapply(names(.rows), function(.e) .ranef[[.e]][.rows[[.e]]])
  names(.eta) <- names(.rows)
  .occ <- list()
  .iov <- .get("iov")
  if (is.list(.iov) && !is.data.frame(.iov)) {
    for (.v in names(.iov)) {
      .d <- .iov[[.v]]
      .rows <- .etaNormalityCols(.d, rows = .etaNormalityMixRows(.ui, .d, .mixNum), drop = .v)
      for (.e in names(.rows)) {
        .x <- .d[[.e]][.rows[[.e]]]
        .lvl <- .d[[.v]][.rows[[.e]]]
        .eta[[.e]] <- .x
        for (.l in unique(.lvl)) {
          .occ[[paste0("rx.", .e, ".", .l)]] <- .x[.lvl == .l]
        }
      }
    }
  }
  list(eta = .eta, occ = .occ)
}

#' Add Shapiro-Wilk rows to a shrinkage table
#'
#' @param shrink a fit's shrinkage table, or a list of them (IOV fits)
#' @param values named list of each column's estimates (a column without
#'   estimates, such as `IWRES`, gets `NA`)
#' @return `shrink` with the `.shrinkShapiroRows` rows appended
#' @author Matthew L. Fidler
#' @noRd
.shrinkAddShapiro <- function(shrink, values) {
  if (is.list(shrink) && !is.data.frame(shrink)) {
    return(lapply(shrink, .shrinkAddShapiro, values = values))
  }
  if (!is.data.frame(shrink) || any(.shrinkShapiroRows %in% row.names(shrink))) {
    return(shrink)
  }
  .add <- as.data.frame(lapply(names(shrink), function(.n) .etaNormalityShapiro(values[[.n]])))
  names(.add) <- names(shrink)
  row.names(.add) <- .shrinkShapiroRows
  rbind(shrink, .add)
}

#' A fit's shrinkage table with its Shapiro-Wilk rows
#'
#' @param fit fit environment
#' @param values the result of `.etaNormalityValues(fit)`
#' @return the shrinkage table (or list of them), or `NULL`
#' @author Matthew L. Fidler
#' @noRd
.shrinkWithShapiro <- function(fit, values = .etaNormalityValues(fit)) {
  .shrink <- get0("shrink", envir = fit, inherits = FALSE)
  if (is.null(.shrink)) {
    return(NULL)
  }
  .shrinkAddShapiro(.shrink, c(values$eta, values$occ))
}

#' One ETA's statistics from a shrinkage table
#'
#' @param shrink a fit's shrinkage table (with Shapiro-Wilk rows), or a list
#'   of them
#' @param eta name of the ETA
#' @return named numeric vector `mean`, `sd`, `meanP`, `shapiroW`, `shapiroP`
#'   (`NA` when absent)
#' @author Matthew L. Fidler
#' @noRd
.etaNormalityShrink <- function(shrink, eta) {
  .rn <- c(
    mean = "mean",
    sd = "sd",
    meanP = "p-value",
    shapiroW = .shrinkShapiroRows[1],
    shapiroP = .shrinkShapiroRows[2]
  )
  .ret <- setNames(rep(NA_real_, length(.rn)), names(.rn))
  .tabs <- if (is.data.frame(shrink)) list(shrink) else shrink
  for (.t in .tabs) {
    if (is.data.frame(.t) && eta %in% names(.t)) {
      .v <- setNames(.t[[eta]], row.names(.t))
      .ret[] <- unname(.v[.rn])
      break
    }
  }
  .ret
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

#' ETA normality table for a fit
#'
#' @param values the result of `.etaNormalityValues()`
#' @param shrink the shrinkage table with its Shapiro-Wilk rows
#' @param skip names of ETAs whose statistics are `NA` (fixed or zero omega)
#' @return data frame with columns `eta`, `n`, `mean`, `sd`, `meanP`,
#'   `shapiroW` and `shapiroP`; one row per ETA
#' @author Matthew L. Fidler
#' @noRd
.etaNormalityDf <- function(values, shrink, skip = character(0)) {
  .eta <- values$eta
  # zero-variance etas are dropped from the fit's etas, so give them NA rows
  .miss <- setdiff(skip, names(.eta))
  .rows <- lapply(c(names(.eta), .miss), function(.e) {
    .s <- .etaNormalityShrink(shrink, .e)
    if (.e %in% skip) {
      .s[] <- NA_real_
    }
    .n <- if (.e %in% .miss) NA_integer_ else sum(is.finite(.eta[[.e]]))
    data.frame(eta = .e, n = .n, as.list(.s), stringsAsFactors = FALSE)
  })
  if (length(.rows) == 0L) {
    return(NULL)
  }
  .ret <- do.call(rbind, .rows)
  rownames(.ret) <- NULL
  .ret
}

#' @rdname nmObjGet
#' @export
nmObjGet.shrink <- function(x, ...) {
  .shrinkWithShapiro(x[[1]])
}
attr(nmObjGet.shrink, "desc") <- "Shrinkage and ETA distribution statistics (incl. Shapiro-Wilk)"

#' @rdname nmObjGet
#' @export
nmObjGet.etaNormality <- function(x, ...) {
  .env <- x[[1]]
  .values <- .etaNormalityValues(.env)
  .ui <- get0("ui", envir = .env, inherits = FALSE)
  if (!is.null(.ui)) {
    .ui <- rxode2::rxUiDecompress(.ui)
  }
  .etaNormalityDf(.values, .shrinkWithShapiro(.env, .values), skip = .etaNormalityFixed(.ui$iniDf))
}
attr(nmObjGet.etaNormality, "desc") <- "ETA zero-mean t-test and Shapiro-Wilk test, from $shrink"
attr(nmObjGet.etaNormality, "rstudio") <- data.frame(
  eta = "eta.cl",
  n = 3L,
  mean = 0,
  sd = 1,
  meanP = 1,
  shapiroW = 1,
  shapiroP = 1
)
