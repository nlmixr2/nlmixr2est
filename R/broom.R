.fixNames <- function(.df) {
  names(.df) <- gsub("Back-.*", "estimate", names(.df))
  names(.df) <- gsub("Estimate", "model.est", names(.df))
  names(.df) <- gsub(".*RSE.*", "rse", names(.df))
  names(.df) <- gsub("SE", "std.error", names(.df))
  names(.df) <- gsub("CI Lower", "conf.low", names(.df))
  names(.df) <- gsub("CI Upper", "conf.high", names(.df))
  names(.df) <- gsub("BSV.*", "bsv", names(.df))
  names(.df) <- gsub("Shrink.*", "shrink", names(.df))
  return(.df)
}

#' A fit's `$parFixedDf` theta rows with broom's column names
#'
#' The table also has rows that are not thetas (a non-mu-referenced eta's
#' variability, an IOV fit's occasion variability, a `bsvcorr()` row); those
#' are not fixed effects and have no estimate to give an interval.
#'
#' @param x nlmixr2 fit
#' @return data.frame of the theta rows, in the table's order; an error when
#'   the table has rows but none is a theta of the model
#' @noRd
.broomParFixedDf <- function(x) {
  .df <- x$parFixedDf
  .ini <- x$ui$iniDf
  .theta <- rownames(.df) %in% .ini$name[!is.na(.ini$ntheta)]
  # every model has a theta, so a table without one is not this fit's
  if (nrow(.df) > 0L && !any(.theta)) {
    stop("no row of $parFixedDf is a theta of the fit's model", call. = FALSE)
  }
  .fixNames(.df[.theta, , drop = FALSE])
}

##' @export
confint.nlmixr2FitCore <- function(object, parm, level = 0.95, ...) {
  .extra <- list(...)
  .exponentiate <- ifelse(any(names(.extra) == "exponentiate"), .extra$exponentiate, FALSE)
  .ciNames <- ifelse(any(names(.extra) == "ciNames"), .extra$ciNames, TRUE)
  .ranPars <- !isFALSE(.extra$ranPars)
  .method <- if (is.null(.extra$method)) "auto" else .extra$method
  .scales <- if (is.null(.extra$scales)) "sdcor" else .extra$scales
  .df <- .broomParFixedDf(object)
  .zv <- qnorm(1 - (1 - level) / 2)
  .low <- .df$model.est - .df$std.error * .zv
  .hi <- .df$model.est + .df$std.error * .zv
  if (is.na(.exponentiate)) {
    .exp <- abs(exp(.df$model.est) - .df$estimate) < 1e-6
    .hi[.exp] <- exp(.hi[.exp])
    .low[.exp] <- exp(.low[.exp])
  } else if (.exponentiate) {
    .hi <- exp(.hi)
    .low <- exp(.low)
  }
  .df <- data.frame(
    model.est = .df$model.est,
    estimate = .df$estimate,
    conf.low = .low,
    conf.high = .hi,
    row.names = rownames(.df)
  )
  # the variability terms follow the thetas, as lme4 lists its variance
  # components; they are on the SD/correlation scale and never exponentiated
  .rp <- if (.ranPars) .ranParTable(object, level, .method, .scales)
  if (!is.null(.rp)) {
    .df <- rbind(
      .df,
      data.frame(
        model.est = .rp$estimate,
        estimate = .rp$estimate,
        conf.low = .rp$conf.low,
        conf.high = .rp$conf.high,
        row.names = .rp$term
      )
    )
  }
  if (!missing(parm)) {
    .df <- .df[parm, , drop = FALSE]
  }
  if (.ciNames) {
    names(.df)[3:4] <- paste(c((1 - level) / 2, (1 - (1 - level) / 2)) * 100, "%")
  }
  .df
}

##' @export
confint.nlmixr2FitCoreSilent <- confint.nlmixr2FitCore

.nlmixr2TidyFixed <- function(x, ..., .ranpar = FALSE) {
  rxode2::rxReq("tibble")
  .extra <- list(...)
  .conf.int <- ifelse(
    any(names(.extra) == "conf.int"),
    .extra$conf.int,
    ifelse(any(names(.extra) == "conf.level"), TRUE, FALSE)
  )
  .conf.level <- ifelse(any(names(.extra) == "conf.level"), .extra$conf.level, 0.95)
  .exponentiate <- ifelse(any(names(.extra) == "exponentiate"), .extra$exponentiate, FALSE)
  .quick <- ifelse(any(names(.extra) == "quick"), .extra$quick, FALSE)
  .rse <- ifelse(any(names(.extra) == "rse"), .extra$rse, FALSE)
  .bsv <- ifelse(any(names(.extra) == "bsv"), .extra$bsv, FALSE)
  .shrink <- ifelse(any(names(.extra) == "shrink"), .extra$shrink, FALSE)
  if (.quick) {
    warning("quick does not do anything for nlmixr2 fit objects")
  }
  .df <- .broomParFixedDf(x)
  .exp <- abs(exp(.df$model.est) - .df$estimate) < 1e-6
  if (is.na(.exponentiate)) {
    ## se(exp(x)) ~= exp(mu_x)*se_x
    ##
    ## Norris N. The Standard Errors of the Geometric and Harmonic Means and Their Application to Index Numbers
    ## Ann. Math. Statist. Volume 11, Number 4 (1940), 445-448.
    .df$std.error[.exp] <- exp(.df$model.est[.exp]) * .df$std.error[.exp]
  } else if (.exponentiate) {
    .df$std.error <- exp(.df$model.est) * .df$std.error
    .df$estimate <- exp(.df$model.est)
  } else if (!.exponentiate) {
    .df$estimate <- .df$model.est
  }
  .df$statistic <- .df$estimate / .df$std.error
  .df$p.value <- stats::pt(.df$statistic, nobs(x) - attr(logLik(x), "df"), lower.tail = FALSE)
  if (.conf.int) {
    .ci <- confint.nlmixr2FitCore(
      x,
      level = .conf.level,
      ciNames = FALSE,
      exponentiate = .exponentiate,
      ranPars = FALSE
    )
    .df$conf.low <- .ci$conf.low
    .df$conf.high <- .ci$conf.high
  } else {
    .df <- .df[, regexpr("^conf[.]", names(.df)) == -1]
  }
  if (!.rse) {
    .df <- .df[, names(.df) != "rse"]
  }
  if (!.bsv) {
    .df <- .df[, names(.df) != "bsv"]
  }
  if (!.shrink) {
    .df <- .df[, names(.df) != "shrink"]
  }
  .ini <- x$ui$iniDf
  # a residual-error theta is a random-effect parameter
  .err <- .ini$err[match(row.names(.df), .ini$name)]
  ## effect   group   term            estimate std.error statistic
  .df <- data.frame(
    effect = "fixed",
    term = row.names(.df),
    .df,
    stringsAsFactors = FALSE
  )
  if (!.ranpar) {
    .df <- .df[is.na(.err), ]
  } else {
    .df <- .df[!is.na(.err), ]
    .df$group <- paste0("Residual(", .err[!is.na(.err)], ")")
    .df$effect <- "ran_pars"
  }
  tibble::as_tibble(.df)
}

.nlmixr2TidyRandom <- function(x, ...) {
  rxode2::rxReq("tibble")
  .omegaR <- x$omegaR
  if (is.list(.omegaR)) {
    .omegaR <- .omegaR$id
  }
  .d <- dim(.omegaR)
  if (.d[1] > 0) {
    .tmp <- stack(x$eta[, -1])
    .df <- data.frame(group = "ID", level = x$eta$ID, term = .tmp$ind, estimate = .tmp$values)
    return(tibble::as_tibble(.df))
  } else {
    return(NULL)
  }
}

#' The fit's `sd__`/`cor__` terms and their Jacobian over the Omega elements
#'
#' One SD per subject-level eta, then one correlation per covariance the model
#' declares (`.omegaDeclaredPairs()`), even one estimated at zero, as
#' `$parFixed`'s `bsvcorr()` rows are.
#'
#' @param iniDf the model's `iniDf`
#' @param omega the subject-level Omega (variance scale)
#' @param omegaR the same Omega as `$omegaR` gives it (SDs and correlations)
#' @param scales `"sdcor"` for `sd__`/`cor__` terms or `"vcov"` for
#'   `var__`/`cov__` terms, as broom.mixed's `scales`
#' @return list of `term`, `estimate`, `isSd` (a diagonal term), `fix` (a
#'   fixed element is involved), `jac` (terms by `om.`/`cov.` elements, the
#'   names `$cov` uses), `alt` (each element's name with its etas swapped),
#'   `sdName`/`elVar` (a variance element's eta name and value, for a `$cov`
#'   that carries it on the SD scale under the eta's own name),
#'   `dep` (which elements each term is a function of) and `varPair` (each
#'   term's two variances, for an off-diagonal term)
#' @noRd
.ranParTerms <- function(iniDf, omega, omegaR, scales = "sdcor") {
  .nm <- dimnames(omega)[[2]]
  .d <- iniDf[!is.na(iniDf$neta1) & iniDf$neta1 == iniDf$neta2, , drop = FALSE]
  .etaNm <- character(max(c(0L, .d$neta1)))
  .etaNm[.d$neta1] <- .d$name
  # in eta-number space; Om only sizes the blocks
  .el <- .omegaDeclaredPairs(diag(length(.etaNm)), iniDf)
  .i <- match(.etaNm[.el[, 1]], .nm)
  .j <- match(.etaNm[.el[, 2]], .nm)
  # an IOV fit's occasion etas are not in the subject-level omega
  .keep <- !is.na(.i) & !is.na(.j)
  .el <- .el[.keep, , drop = FALSE]
  .i <- .i[.keep]
  .j <- .j[.keep]
  .elFix <- .omegaFixed(iniDf, .el)
  .diag <- .i == .j
  .sdE <- match(seq_along(.nm), ifelse(.diag, .i, NA_integer_))
  .offE <- which(!.diag)
  .offE <- .offE[order(pmin(.i[.offE], .j[.offE]), pmax(.i[.offE], .j[.offE]))]
  .a <- pmax(.i[.offE], .j[.offE])
  .b <- pmin(.i[.offE], .j[.offE])
  .sep1 <- getOption("broom.mixed.sep1", "__")
  .sep2 <- getOption("broom.mixed.sep2", ".")
  .pre <- if (scales == "vcov") c("var", "cov") else c("sd", "cor")
  .term <- c(
    paste0(.pre[1], .sep1, .nm),
    if (length(.offE) == 0) character(0) else paste0(.pre[2], .sep1, .nm[.a], .sep2, .nm[.b])
  )
  .nSd <- length(.nm)
  .jac <- matrix(0, length(.term), nrow(.el), dimnames = list(.term, .foceiOmegaCovNames(.el, .etaNm)))
  .dep <- .jac != 0
  .sd <- sqrt(diag(omega))
  .varPair <- rbind(cbind(diag(omega), diag(omega)), cbind(diag(omega)[.a], diag(omega)[.b]))
  .alt <- .foceiOmegaCovNames(.el[, 2:1, drop = FALSE], .etaNm)
  .sdName <- ifelse(.diag, .etaNm[.el[, 1]], NA_character_)
  .elVar <- ifelse(.diag, diag(omega)[.i], NA_real_)
  if (scales == "vcov") {
    # the terms are the elements themselves
    .e <- c(.sdE, .offE)
    .jac[cbind(seq_along(.e), .e)] <- 1
    return(list(
      term = .term,
      estimate = c(diag(omega), omega[cbind(.a, .b)]),
      isSd = seq_along(.term) <= .nSd,
      fix = .elFix[.e],
      jac = .jac,
      dep = .jac != 0,
      alt = .alt,
      sdName = .sdName,
      elVar = .elVar,
      varPair = .varPair
    ))
  }
  for (.k in seq_len(.nSd)) {
    .jac[.k, .sdE[.k]] <- 1 / (2 * .sd[.k])
    .dep[.k, .sdE[.k]] <- TRUE
  }
  for (.m in seq_along(.offE)) {
    .k <- .nSd + .m
    .r <- omega[.a[.m], .b[.m]] / (.sd[.a[.m]] * .sd[.b[.m]])
    .e <- c(.offE[.m], .sdE[.a[.m]], .sdE[.b[.m]])
    .jac[.k, .e] <- c(
      1 / (.sd[.a[.m]] * .sd[.b[.m]]),
      -.r / (2 * omega[.a[.m], .a[.m]]),
      -.r / (2 * omega[.b[.m], .b[.m]])
    )
    .dep[.k, .e] <- TRUE
  }
  list(
    term = .term,
    estimate = c(diag(omegaR), omegaR[cbind(.a, .b)]),
    isSd = seq_along(.term) <= .nSd,
    fix = as.vector(.dep %*% .elFix > 0),
    jac = .jac,
    dep = .dep,
    alt = .alt,
    sdName = .sdName,
    elVar = .elVar,
    varPair = .varPair
  )
}

#' SEs of the `sd__`/`cor__` terms from a fit's Omega covariance
#'
#' @param x nlmixr2 fit
#' @param terms `.ranParTerms()`
#' @return SEs, the square root of the diagonal of `J cov J'`; `NA` where
#'   `x$cov` lacks an element the term is a function of
#' @noRd
.ranParSe <- function(x, terms) {
  .cov <- x$cov
  .rn <- if (is.matrix(.cov)) rownames(.cov)
  .cn <- colnames(terms$jac)
  # each element as $cov names it: om./cov. in either eta order, or (FOCEi's
  # IOV magnitude) the eta's own name on the SD scale, d var / d sd = 2 sd
  .use <- ifelse(.cn %in% .rn, .cn, ifelse(terms$alt %in% .rn, terms$alt, NA_character_))
  .sdRow <- is.na(.use) & terms$sdName %in% .rn
  .use[.sdRow] <- terms$sdName[.sdRow]
  .fac <- ifelse(.sdRow, 2 * sqrt(terms$elVar), 1)
  # one term at a time over its own elements, so an NA row of a partial
  # (rank-deficient) $cov reaches only the terms that use it
  .v <- vapply(
    seq_len(nrow(terms$jac)),
    function(k) {
      .d <- terms$dep[k, ]
      if (!any(.d) || anyNA(.use[.d])) {
        return(NA_real_)
      }
      .c <- .use[.d]
      .covJacobian(.cov[.c, .c, drop = FALSE], terms$jac[k, .d, drop = FALSE] * .fac[.d], NULL)[1, 1]
    },
    numeric(1)
  )
  ifelse(is.finite(.v) & .v >= 0, sqrt(pmax(.v, 0)), NA_real_)
}

#' Wald intervals of the variability terms
#'
#' SDs and variances on the log scale and correlations on atanh, as
#' nlme::intervals(), so each stays in range; covariances on their own scale.
#'
#' @param pars term estimates
#' @param se their SEs
#' @param isSd which terms are diagonal (SDs or variances)
#' @param level confidence level
#' @param scales `"sdcor"` or `"vcov"`
#' @return list of `low` and `high`
#' @noRd
.ranParCi <- function(pars, se, isSd, level, scales = "sdcor") {
  .zse <- qnorm(1 - (1 - level) / 2) * se
  .low <- .high <- rep(NA_real_, length(pars))
  .s <- pars[isSd]
  .low[isSd] <- .s * exp(-.zse[isSd] / .s)
  .high[isSd] <- .s * exp(.zse[isSd] / .s)
  .r <- pars[!isSd]
  if (scales == "vcov") {
    .low[!isSd] <- .r - .zse[!isSd]
    .high[!isSd] <- .r + .zse[!isSd]
  } else {
    .d <- .zse[!isSd] / (1 - .r^2)
    .low[!isSd] <- tanh(atanh(.r) - .d)
    .high[!isSd] <- tanh(atanh(.r) + .d)
  }
  list(low = .low, high = .high)
}

#' SEs and intervals of the variability terms from the Wishart distribution
#'
#' Omega estimated from `n` subjects' etas is taken as Wishart with `n`
#' degrees of freedom: each variance is a scaled chi-square, each correlation
#' uses Fisher's z with SE `1/sqrt(n - 3)` and each covariance `w_ab` has
#' variance `(w_ab^2 + w_aa w_bb) / n`.
#'
#' @inheritParams .ranParCi
#' @param n number of subjects
#' @param varPair each term's two variances (`.ranParTerms()`), for `"vcov"`
#' @return list of `se`, `low` and `high`
#' @noRd
.ranParWishart <- function(pars, isSd, n, level, scales = "sdcor", varPair = NULL) {
  .a <- (1 - level) / 2
  .se <- .low <- .high <- rep(NA_real_, length(pars))
  .s <- pars[isSd]
  # the chi-square interval of a variance, on the term's own scale
  .pw <- if (scales == "vcov") 1 else 0.5
  .se[isSd] <- if (scales == "vcov") .s * sqrt(2 / n) else .s / sqrt(2 * n)
  .low[isSd] <- .s * (n / stats::qchisq(1 - .a, n))^.pw
  .high[isSd] <- .s * (n / stats::qchisq(.a, n))^.pw
  .r <- pars[!isSd]
  if (scales == "vcov") {
    .vp <- varPair[!isSd, , drop = FALSE]
    .se[!isSd] <- sqrt((.r^2 + .vp[, 1] * .vp[, 2]) / n)
    .low[!isSd] <- .r - stats::qnorm(1 - .a) * .se[!isSd]
    .high[!isSd] <- .r + stats::qnorm(1 - .a) * .se[!isSd]
  } else {
    .se[!isSd] <- (1 - .r^2) / sqrt(n)
    if (isTRUE(n > 3)) {
      .d <- stats::qnorm(1 - .a) / sqrt(n - 3)
      .low[!isSd] <- tanh(atanh(.r) - .d)
      .high[!isSd] <- tanh(atanh(.r) + .d)
    }
  }
  list(se = .se, low = .low, high = .high)
}

#' The fit's `sd__`/`cor__` terms with SEs and intervals
#'
#' @param x nlmixr2 fit
#' @param level confidence level
#' @param method `"auto"` (the Omega rows of `$cov` when present, else
#'   Wishart), `"cov"` or `"wishart"`
#' @param scales `"sdcor"` or `"vcov"` (`.ranParTerms()`)
#' @return data.frame of `group` (`"ID"`, or an IOV level's occasion
#'   variable), `term`, `estimate`, `std.error`, `conf.low`, `conf.high`,
#'   `method`; `NULL` when the model has no etas
#' @noRd
.ranParTable <- function(x, level = 0.95, method = c("auto", "cov", "wishart"), scales = c("sdcor", "vcov")) {
  method <- match.arg(method)
  scales <- match.arg(scales)
  .omega <- x$omega
  .omegaR <- x$omegaR
  if (!is.list(.omega)) {
    .omega <- list(id = .omega)
    .omegaR <- list(id = .omegaR)
  }
  # the model as estimated, as $parFixed reports its correlations: a ui
  # rebuilt from the estimated omega has lost a covariance estimated at zero
  .ini0 <- if (is.environment(x$env)) .updateParFixedBsvCorIni(x$env, x$ui) else x$ui$iniDf
  .ret <- lapply(names(.omega), function(.lvl) {
    .om <- .omega[[.lvl]]
    if (is.null(.om) || nrow(.om) == 0) {
      return(NULL)
    }
    # an IOV fit's estimated model names its occasion etas rx.<eta>; the
    # final ui names them as $omega does
    .ini <- if (all(dimnames(.om)[[2]] %in% .ini0$name[!is.na(.ini0$neta1)])) .ini0 else x$ui$iniDf
    # the number of subjects, or of subject-occasions for an IOV level
    .n <- if (.lvl == "id") {
      nrow(x$eta)
    } else if (is.environment(x$env) && is.data.frame(x$env$iov[[.lvl]])) {
      nrow(x$env$iov[[.lvl]])
    } else {
      NA_integer_
    }
    .df <- .ranParLevel(x, .ini, .om, .omegaR[[.lvl]], .n, level, method, scales)
    data.frame(group = if (.lvl == "id") "ID" else .lvl, .df, stringsAsFactors = FALSE)
  })
  .ret <- do.call(rbind, .ret)
  if (is.null(.ret) || nrow(.ret) == 0) {
    return(NULL)
  }
  rownames(.ret) <- NULL
  .ret
}

#' One Omega level's terms for `.ranParTable()`
#'
#' @inheritParams .ranParTable
#' @param ini the model's `iniDf`
#' @param omega,omegaR the level's Omega and its SD/correlation form
#' @param n subjects (or subject-occasions) for the Wishart intervals
#' @return data.frame as `.ranParTable()`, without `group`
#' @noRd
.ranParLevel <- function(x, ini, omega, omegaR, n, level, method, scales) {
  .t <- .ranParTerms(ini, omega, omegaR, scales)
  .pars <- setNames(.t$estimate, .t$term)
  .isSd <- .t$isSd
  .cse <- if (method == "wishart") rep(NA_real_, length(.pars)) else .ranParSe(x, .t)
  .ci <- .ranParCi(.pars, .cse, .isSd, level, scales)
  .useW <- if (method == "cov") rep(FALSE, length(.pars)) else is.na(.cse)
  .w <- .ranParWishart(.pars, .isSd, n, level, scales, .t$varPair)
  .df <- data.frame(
    term = names(.pars),
    estimate = unname(.pars),
    std.error = ifelse(.useW, .w$se, .cse),
    conf.low = ifelse(.useW, .w$low, .ci$low),
    conf.high = ifelse(.useW, .w$high, .ci$high),
    method = ifelse(.useW, "wishart", "cov"),
    stringsAsFactors = FALSE
  )
  # a fixed Omega element is not estimated, nor is any under est="posthoc"
  .fix <- .t$fix | identical(x$est, "posthoc")
  .df[.fix, c("std.error", "conf.low", "conf.high")] <- NA_real_
  .df$method[.fix | is.na(.df$std.error)] <- NA_character_
  .df
}

.nlmixr2TidyRandomPar <- function(x, ...) {
  rxode2::rxReq("tibble")
  .extra <- list(...)
  .level <- if (is.null(.extra$conf.level)) 0.95 else .extra$conf.level
  .method <- if (is.null(.extra$conf.method)) "auto" else .extra$conf.method
  .scales <- if (is.null(.extra$.ranScale)) "sdcor" else .extra$.ranScale
  .rp <- .ranParTable(x, .level, .method, .scales)
  if (!is.null(.rp)) {
    .p1 <- data.frame(
      effect = "ran_pars",
      group = .rp$group,
      term = .rp$term,
      estimate = .rp$estimate,
      std.error = .rp$std.error,
      statistic = NA_real_,
      p.value = NA_real_,
      stringsAsFactors = FALSE
    )
    if (isTRUE(.extra$conf.int) || (is.null(.extra$conf.int) && !is.null(.extra$conf.level))) {
      .p1$conf.low <- .rp$conf.low
      .p1$conf.high <- .rp$conf.high
    }
    .p1 <- .reorderCols(.p1)
    # broom.mixed exponentiates fixed effects only
    .args <- c(list(x), .extra[names(.extra) != "exponentiate"], list(.ranpar = TRUE))
    .p2 <- data.frame(do.call(.nlmixr2TidyFixed, .args), stringsAsFactors = FALSE) |>
      .reorderCols()
    .df <- rbind(.p1, .p2)
    for (.v in c("statistic", "p.value", "std.error")) {
      if (all(is.na(.df[, .v]))) {
        .df <- .df[, names(.df) != .v]
      }
    }
    return(tibble::as_tibble(.df))
  } else {
    return(NULL)
  }
}
##   effect   group    term                  estimate std.error statistic
##   <chr>    <chr>    <chr>                    <dbl>     <dbl>     <dbl>
## 1 fixed    NA       (Intercept)           251.          6.82     36.8
## 2 fixed    NA       Days                   10.5         1.55      6.77
## 3 ran_pars Subject  sd__(Intercept)        24.7        NA        NA
## 4 ran_pars Subject  sd__Days                5.92       NA        NA
## 5 ran_pars Subject  cor__(Intercept).Days   0.0656     NA        NA
## 6 ran_pars Residual sd__Observation        25.6        NA        NA

## Row names and order taken & adapted from
## https://github.com/bbolker/broom.mixed/blob/master/R/utilities.R#L238-L248
.reorderCols <- function(x) {
  allCols <- c(
    "response",
    "effect",
    "component", ## glmmTMB, brms
    "group",
    "level",
    "term",
    "index",
    "estimate",
    "std.error",
    "statistic",
    "df",
    "p.value",
    "conf.low",
    "conf.high",
    "rhat",
    "ess"
  )
  return(x[, intersect(allCols, names(x))])
}

.coefPar <- function(x, exponentiate = FALSE, ...) {
  rxode2::rxReq("tibble")
  .omegaR <- x$omegaR
  if (is.list(.omegaR)) {
    .omegaR <- .omegaR$id
  }
  .d <- dim(.omegaR)
  if (.d[1] == 0) {
    return(NULL)
  }
  .muRef <- x$ui$muRefDataFrame
  .curEval <- x$ui$muRefCurEval
  .theta <- x$theta
  .df <- .broomParFixedDf(x)
  if (is.na(exponentiate)) {
    .exp <- abs(exp(.df$model.est) - .df$estimate) < 1e-6
  } else if (exponentiate) {
    .exp <- setNames(rep(TRUE, length(.theta)), names(.theta))
  } else {
    .exp <- setNames(rep(FALSE, length(.theta)), names(.theta))
  }
  .eta <- x$eta
  .x <- setNames(
    data.frame(lapply(names(.eta), function(eta) {
      .w <- which(.muRef$eta == eta)
      if (length(.w) == 1L) {
        .thetaName <- .muRef$theta[.w]
        .ret <- .eta[[eta]] + .theta[.thetaName]
        if (any(.thetaName == names(.exp))) {
          if (.exp[.thetaName]) {
            .ret <- exp(.ret)
          }
        }
      } else {
        if (eta != "ID") {
          warning(
            sprintf("the parameter '%s' is not mu-referenced and the coef will not be returned", eta),
            call. = FALSE
          )
        }
        .ret <- .eta[[eta]]
      }
      return(.ret)
    })),
    sapply(names(.eta), function(eta) {
      .w <- which(.muRef$eta == eta)
      if (length(.w) == 1L) {
        .thetaName <- .muRef$theta[.w]
        return(.thetaName)
      }
      return(eta)
    })
  )
  .tmp <- stack(.x[, -1])
  .df <- data.frame(group = "ID", level = .x$ID, term = .tmp$ind, estimate = .tmp$values)
  return(tibble::as_tibble(.df))
}

tidy.nlmixr2FitCore <- function(x, ...) {
  rxode2::rxReq("tibble")
  rxode2::rxReq("dplyr")
  .extra <- list(...)
  if (any(names(.extra) == "effects")) {
    .effects <- .extra$effects
  } else {
    if (any(names(.extra) == "effect")) {
      .effects <- .extra$effect
    } else {
      .effects <- c("fixed", "ran_pars")
    }
  }
  .effects <- match.arg(
    .effects,
    c("fixed", "random", "ran_vals", "ran_pars", "ran_coef", "ran_coefs"),
    several.ok = TRUE
  )
  # one scale per effect, as broom.mixed; only ran_pars has one
  .scales <- .extra$scales
  .ranScale <- "sdcor"
  if (!is.null(.scales)) {
    if (length(.scales) != length(.effects)) {
      stop("'scales' needs one value (or NA) per element of 'effects'", call. = FALSE)
    }
    .s <- .scales[.effects == "ran_pars"]
    if (length(.s) == 1L && !is.na(.s)) {
      .ranScale <- match.arg(.s, c("sdcor", "vcov"))
    }
  }
  .ret <- list()
  if (any(.effects == "fixed")) {
    .ret$fixed <- .nlmixr2TidyFixed(x, ...)
  }
  if (any(.effects == "random") || any(.effects == "ran_vals")) {
    .ret$ran_vals <- .nlmixr2TidyRandom(x, ...)
  }
  if (any(.effects == "ran_pars")) {
    .ret$ran_pars <- .nlmixr2TidyRandomPar(x, ..., .ranScale = .ranScale)
  }
  # broom.mixed spells it ran_coefs; ran_coef is kept for older code
  for (.e in intersect(c("ran_coef", "ran_coefs"), .effects)) {
    .ret[[.e]] <- .coefPar(x, ...)
  }
  if (all(unlist(lapply(seq_along(.ret), is.null)))) {
    return(NULL)
  }
  dplyr::bind_rows(.ret, .id = "effect") |>
    tibble::as_tibble() |>
    .reorderCols()
}

tidy.nlmixr2FitCoreSilent <- tidy.nlmixr2FitCore

glance.nlmixr2FitCore <- function(x, ...) {
  rxode2::rxReq("tibble")
  .lst <- list(...)
  if (any(names(.lst) == "type")) {
    setOfv(x, type = .lst$type)
  }
  .aic <- AIC(x) ## To calculate AIC if needed
  .df <- x$objDf
  if (length(x$logLik) > 1L) {
    .df <- .df[which(tolower(rownames(x$objDf)) == tolower(x$ofvType)), ]
  }
  names(.df) <- gsub("Log-likelihood", "logLik", names(.df))
  names(.df) <- gsub("Condition#(Cor)", "conditionNumberCor", names(.df))
  names(.df) <- gsub("Condition#(Cov)", "conditionNumberCov", names(.df))
  tibble::as_tibble(.df)
}

glance.nlmixr2FitCoreSilent <- glance.nlmixr2FitCore

augment.nlmixr2FitCore <- function(x, ...) {
  stop("augment is not yet implemented for nlmixr2 models", call. = FALSE)
}

augment.nlmixr2FitCoreSilent <- augment.nlmixr2FitCore
