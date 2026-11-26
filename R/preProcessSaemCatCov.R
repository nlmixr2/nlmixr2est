# Pre-processing hook: mu-reference multiplicative 0/1 covariate effects for saem
#
# `cl <- exp(tcl + eta.cl) * (1 + clSex * SEX)` puts `clSex` outside the
# mu-referenced phi, where saem's conditional M-step cannot separate it from
# `eta.cl` (issue #1315).  For a covariate that only takes the values 0 and 1,
# `1 + clSex*SEX == exp(log(1 + clSex)*SEX)` exactly, so the line is fit as
# `exp(tcl + eta.cl + rxBoundedTr.clSex*SEX)` and `clSex` is back-transformed
# with the bounded-transform machinery (type "log1p").

#' Is a data column a 0/1 covariate?
#'
#' @param data estimation data
#' @param cov column name
#' @return `TRUE` when every non-missing value is 0 or 1
#' @noRd
.saemCatCovIsBinary <- function(data, cov) {
  if (!is.data.frame(data) || !(cov %in% names(data))) {
    return(FALSE)
  }
  .v <- data[[cov]]
  if (!is.numeric(.v) && !is.logical(.v)) {
    return(FALSE)
  }
  .v <- .v[!is.na(.v)]
  length(.v) > 0L && all(.v %in% c(0, 1))
}

#' Flatten a product into its factors
#'
#' @param e language object
#' @return list of factors (a bare product inside parentheses is flattened too)
#' @noRd
.saemCatCovFactors <- function(e) {
  if (is.call(e) && identical(e[[1]], quote(`(`))) {
    .in <- e[[2]]
    if (is.call(.in) && identical(.in[[1]], quote(`*`))) {
      return(.saemCatCovFactors(.in))
    }
    return(list(e))
  }
  if (is.call(e) && identical(e[[1]], quote(`*`)) && length(e) == 3L) {
    return(c(.saemCatCovFactors(e[[2]]), .saemCatCovFactors(e[[3]])))
  }
  list(e)
}

#' Match a `(1 + theta*COV)` factor
#'
#' @param f factor
#' @param thetas candidate theta names
#' @param covs 0/1 covariate names
#' @return `c(theta=, cov=)` or `NULL`
#' @noRd
.saemCatCovMatch <- function(f, thetas, covs) {
  while (is.call(f) && identical(f[[1]], quote(`(`))) {
    f <- f[[2]]
  }
  if (!is.call(f) || !identical(f[[1]], quote(`+`)) || length(f) != 3L) {
    return(NULL)
  }
  .one <- function(x) is.numeric(x) && length(x) == 1L && x == 1
  .term <- if (.one(f[[2]])) f[[3]] else if (.one(f[[3]])) f[[2]] else NULL
  if (!is.call(.term) || !identical(.term[[1]], quote(`*`)) || length(.term) != 3L) {
    return(NULL)
  }
  .a <- .term[[2]]
  .b <- .term[[3]]
  if (!is.name(.a) || !is.name(.b)) {
    return(NULL)
  }
  .a <- as.character(.a)
  .b <- as.character(.b)
  if (.a %in% thetas && .b %in% covs) {
    return(c(theta = .a, cov = .b))
  }
  if (.b %in% thetas && .a %in% covs) {
    return(c(theta = .b, cov = .a))
  }
  NULL
}

#' Rewrite one model line's `(1 + theta*COV)` factors into its `exp()` factor
#'
#' @param expr model line
#' @param thetas candidate theta names
#' @param covs 0/1 covariate names
#' @param etas eta names
#' @return list of `expr` (rewritten line) and `match` (matrix of theta/cov),
#'   or `NULL` when the line has no such factor
#' @noRd
.saemCatCovLine <- function(expr, thetas, covs, etas) {
  if (!is.call(expr) || !(identical(expr[[1]], quote(`<-`)) || identical(expr[[1]], quote(`=`)))) {
    return(NULL)
  }
  .fac <- .saemCatCovFactors(expr[[3]])
  if (length(.fac) < 2L) {
    return(NULL)
  }
  .isExp <- vapply(
    .fac,
    function(f) {
      is.call(f) && identical(f[[1]], quote(exp)) && length(f) == 2L && any(all.vars(f) %in% etas)
    },
    logical(1)
  )
  if (sum(.isExp) != 1L) {
    return(NULL)
  }
  .m <- lapply(.fac, .saemCatCovMatch, thetas = thetas, covs = covs)
  .cat <- !vapply(.m, is.null, logical(1))
  if (!any(.cat)) {
    return(NULL)
  }
  .arg <- .fac[[which(.isExp)]][[2]]
  for (.k in which(.cat)) {
    .arg <- bquote(.(.arg) + .(as.name(paste0("rxBoundedTr.", .m[[.k]][["theta"]]))) * .(as.name(.m[[.k]][["cov"]])))
  }
  .fac[[which(.isExp)]] <- bquote(exp(.(.arg)))
  .fac <- .fac[!.cat]
  .rhs <- Reduce(function(a, b) bquote(.(a) * .(b)), .fac)
  .new <- expr
  .new[[3]] <- .rhs
  list(expr = .new, match = do.call(rbind, .m[.cat]))
}

#' Number of times each name appears in the model
#'
#' @param lstExpr model lines
#' @param names names to count
#' @return named integer vector
#' @noRd
.saemCatCovCount <- function(lstExpr, names) {
  .all <- unlist(lapply(lstExpr, function(e) all.names(e)))
  vapply(names, function(n) sum(.all == n), integer(1))
}

#' Build the saem ui with the selected factors mu-referenced
#'
#' @param ui rxode2 ui
#' @param rew list of line rewrites (line index, rewrite)
#' @param thetas thetas to rewrite
#' @return rewritten ui carrying `boundedTransforms`
#' @noRd
.saemCatCovUi <- function(ui, rew, thetas) {
  .lstExpr <- ui$lstExpr
  .iniDf <- ui$iniDf
  .transforms <- list()
  for (.r in rew) {
    .keep <- .r$match[, "theta"] %in% thetas
    if (!any(.keep)) {
      next
    }
    .line <- .saemCatCovLine(
      .lstExpr[[.r$i]],
      .r$match[.keep, "theta"],
      .r$match[.keep, "cov"],
      .r$etas
    )
    for (.t in .line$match[, "theta"]) {
      .w <- which(.iniDf$name == .t)
      .transforms[[length(.transforms) + 1L]] <- list(
        name = .t,
        internalName = paste0("rxBoundedTr.", .t),
        type = "log1p",
        lower = .iniDf$lower[.w],
        upper = .iniDf$upper[.w],
        initTrans = log1p(.iniDf$est[.w]),
        initOrig = .iniDf$est[.w],
        lineOrig = .lstExpr[[.r$i]],
        lineNew = .line$expr
      )
      .iniDf[.w, c("name", "lower", "upper", "est")] <-
        list(paste0("rxBoundedTr.", .t), -Inf, Inf, log1p(.iniDf$est[.w]))
    }
    .lstExpr[[.r$i]] <- .line$expr
  }
  .model <- str2lang(paste0("model({", paste(vapply(.lstExpr, deparse1, character(1)), collapse = "\n"), "})"))
  .ini <- as.expression(lotri::as.lotri(.iniDf))
  .ini[[1]] <- quote(`ini`)
  .fun <- .getUiFunFromIniAndModel(ui, .ini, .model)
  .newUi <- rxode2::rxUiDecompress(suppressWarnings(.fun()))
  assign("modelName", ui$modelName, envir = .newUi)
  .newUi$boundedTransforms <- c(ui$boundedTransforms, .transforms)
  .newUi
}

#' Pre-processing hook: fit `exp(mu + eta) * (1 + theta*COV)` mu-referenced
#'
#' Only for `est="saem"` with `muRefCovAlg=TRUE`, a 0/1 data covariate, a
#' free structural theta used nowhere else and a start above -1.
#'
#' @param ui rxode2 ui
#' @param est estimation method
#' @param data estimation data
#' @param control saem control
#' @return list with the rewritten ui, or `NULL` when nothing changes
#' @noRd
.preProcessSaemCatCov <- function(ui, est, data, control) {
  if (!identical(est, "saem") || !isTRUE(control$muRefCovAlg)) {
    return(NULL)
  }
  ui <- rxode2::rxUiDecompress(ui)
  if (length(tryCatch(ui$nlmixr2LikPenalty, error = function(e) NULL)) > 0L) {
    return(NULL)
  }
  .iniDf <- ui$iniDf
  .covs <- ui$allCovs
  .covs <- .covs[vapply(.covs, .saemCatCovIsBinary, logical(1), data = data)]
  if (length(.covs) == 0L) {
    return(NULL)
  }
  .th <- .iniDf[!is.na(.iniDf$ntheta) & !.iniDf$fix & is.na(.iniDf$err) & .iniDf$est > -1, , drop = FALSE]
  .thetas <- .th$name[!grepl("^rxBoundedTr[.]", .th$name)]
  .etas <- .iniDf$name[!is.na(.iniDf$neta1) & .iniDf$neta1 == .iniDf$neta2]
  if (length(.thetas) == 0L || length(.etas) == 0L) {
    return(NULL)
  }
  .lstExpr <- ui$lstExpr
  .rew <- list()
  for (.i in seq_along(.lstExpr)) {
    .r <- .saemCatCovLine(.lstExpr[[.i]], .thetas, .covs, .etas)
    if (!is.null(.r)) {
      .rew[[length(.rew) + 1L]] <- list(i = .i, match = .r$match, etas = .etas)
    }
  }
  if (length(.rew) == 0L) {
    return(NULL)
  }
  .cand <- unlist(lapply(.rew, function(r) r$match[, "theta"]))
  # a theta used anywhere else cannot move into the mu-referenced phi
  .cand <- .cand[.saemCatCovCount(.lstExpr, .cand) == 1L]
  if (length(.cand) == 0L) {
    return(NULL)
  }
  .newUi <- .saemCatCovUi(ui, .rew, .cand)
  .ok <- .cand[paste0("rxBoundedTr.", .cand) %in% .newUi$muRefCovariateDataFrame$covariateParameter]
  if (length(.ok) == 0L) {
    return(NULL)
  }
  if (length(.ok) < length(.cand)) {
    .newUi <- .saemCatCovUi(ui, .rew, .ok)
  }
  for (.tr in .newUi$boundedTransforms) {
    if (!identical(.tr$type, "log1p")) {
      next
    }
    warning(sprintf("'%s' fit as log(1+%s), mu-referenced", .tr$name, .tr$name), call. = FALSE)
    if (.tr$lower > -1 || is.finite(.tr$upper)) {
      warning(sprintf("bounds of '%s' not enforced by saem", .tr$name), call. = FALSE)
    }
  }
  list(ui = .newUi)
}

preProcessHooksAdd(".preProcessSaemCatCov", .preProcessSaemCatCov)
