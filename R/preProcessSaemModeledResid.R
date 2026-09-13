#' Normal endpoints with a modeled (non-theta) residual error component
#'
#' rxode2 stores a model variable used as an error argument (`a` in
#' `a <- add.sd*exp(eta.sd); cp ~ add(a)`) by name in the `predDf` argument
#' columns; a plain theta leaves them `NA`.
#'
#' @param ui rxode2 ui
#' @return character vector of endpoint conditions
#' @noRd
.saemModeledResidualCond <- function(ui) {
  .pred <- ui$predDf
  if (is.null(.pred) || length(.pred$cond) == 0L) return(character(0))
  .cols <- intersect(c("a", "b", "c", "d", "e", "f", "lambda"), names(.pred))
  if (length(.cols) == 0L) return(character(0))
  .modeled <- vapply(seq_along(.pred$cond), function(i) {
    .pred$distribution[i] == "norm" &&
      any(!is.na(unlist(.pred[i, .cols, drop = TRUE])))
  }, logical(1), USE.NAMES = FALSE)
  as.character(.pred$cond[.modeled])
}

#' Append `+ dnorm()` to an error line, before any `| condition`
#'
#' @param line error model line, e.g. `cp ~ add(a) | cp`
#' @return the line with `dnorm()` added
#' @noRd
.saemAddDnormToErrLine <- function(line) {
  .rhs <- line[[3]]
  if (is.call(.rhs) && identical(.rhs[[1]], as.name("|"))) {
    .rhs[[2]] <- call("+", .rhs[[2]], quote(dnorm()))
  } else {
    .rhs <- call("+", .rhs, quote(dnorm()))
  }
  line[[3]] <- .rhs
  line
}

#' Valid range of a likelihood argument held in a `predDf` column
#'
#' Uses rxode2's own argument ranges (the ones it applies to a theta passed
#' directly as that argument).
#'
#' @param pred1 one `predDf` row
#' @param col the `predDf` column holding the argument
#' @return list with `range` (`c(lower, upper)` or `NULL` when unknown) and
#'   `skip` (`TRUE` for a count argument such as the binomial size)
#' @noRd
.saemPredArgRange <- function(pred1, col) {
  .ns <- asNamespace("rxode2")
  .map <- get(".namedArgumentsToPredDf", envir = .ns)
  .ranges <- get(".errDistArgRanges", envir = .ns)
  .dist <- paste(pred1$distribution)
  .keys <- if (.dist %in% c("norm", "dnorm", "t", "cauchy")) {
    c(strsplit(paste(pred1$errType), " + ", fixed = TRUE)[[1]], .dist)
  } else {
    .dist
  }
  for (.key in .keys) {
    .pos <- match(col, .map[[.key]])
    if (is.na(.pos)) next
    .tag <- if (.pos == 1L) .key else paste0(.key, .pos)
    if (.tag == "binom") return(list(range = NULL, skip = TRUE))
    return(list(range = .ranges[[.tag]], skip = FALSE))
  }
  list(range = NULL, skip = FALSE)
}

#' The theta a model variable is a plain alias of
#'
#' Follows `a <- b; b <- theta` style assignments that precede line `idx`.
#'
#' @param lines deparsed model lines
#' @param lhs `.rxLineLhs(lines)`
#' @param idx line index of the error model
#' @param v variable name
#' @return the aliased name, or `NA` when `v` is an expression of several terms
#' @noRd
.saemAliasTheta <- function(lines, lhs, idx, v) {
  .cur <- v
  for (.k in seq_len(50L)) {
    .w <- which(!is.na(lhs) & lhs == .cur & seq_along(lhs) < idx)
    if (length(.w) == 0L) return(.cur)
    .rhs <- trimws(sub(";[ \t]*$", "", sub("^[^<=~]*(<-|=|~)", "", lines[max(.w)])))
    if (!grepl("^[A-Za-z._][A-Za-z0-9._]*$", .rhs)) return(NA_character_)
    .cur <- .rhs
  }
  NA_character_
}

#' Estimated likelihood thetas without an eta in a saem general-likelihood fit
#'
#' saem refines an eta-less theta only through its phi0 step, which leaves a
#' likelihood parameter near its initial value; these thetas get a temporary
#' eta.  A theta that is the likelihood argument itself uses that argument's
#' range; any other theta uses its own `ini()` bounds.  Thetas that feed a
#' prediction, count arguments, transform lambdas and `ll()`/ordinal endpoints
#' are left alone.
#'
#' @param ui rxode2 ui
#' @return data frame with `theta`, `lower` and `upper`
#' @noRd
.saemPseudoEtaThetas <- function(ui) {
  .empty <- data.frame(theta = character(0), lower = numeric(0), upper = numeric(0))
  if (!.saemGeneralLik(ui)) return(.empty)
  .pred <- ui$predDf
  .iniDf <- ui$iniDf
  .lines <- vapply(ui$lstExpr, deparse1, character(1))
  .lhs <- .rxLineLhs(.lines)
  .cols <- intersect(c("a", "b", "c", "d", "e", "f"), names(.pred))
  .resid <- character(0)
  .predDeps <- character(0)
  .skip <- character(0)
  .llThetas <- character(0)
  .argRange <- list()
  for (i in seq_along(.pred$cond)) {
    if (paste(.pred$distribution[i]) == "ordinal") next
    .idx <- .pred$line[i]
    if (paste(.pred$distribution[i]) == "LL") {
      # ll() has no separate prediction, and a linCmt()/ODE state hides the
      # assignments behind it, so every theta in the model informs it
      .llThetas <- unique(c(.llThetas, unlist(lapply(ui$lstExpr, all.vars))))
      next
    }
    .predDeps <- c(.predDeps, .rxMtimeDeps(.lines, .lhs, .idx, paste(.pred$var[i])))
    for (.c in .cols) {
      .v <- .pred[[.c]][i]
      if (is.na(.v)) next
      .resid <- c(.resid, .rxMtimeDeps(.lines, .lhs, .idx, .v))
      .alias <- .saemAliasTheta(.lines, .lhs, .idx, .v)
      if (is.na(.alias)) next
      .r <- .saemPredArgRange(.pred[i, , drop = FALSE], .c)
      if (.r$skip) .skip <- c(.skip, .alias)
      else if (!is.null(.r$range)) .argRange[[.alias]] <- .r$range
    }
    for (.w in which(!is.na(.iniDf$err) & .iniDf$condition == .pred$cond[i])) {
      if (grepl("boxCox|yeoJohnson|^ar$|^binom$", .iniDf$err[.w])) {
        .skip <- c(.skip, .iniDf$name[.w])
      } else {
        .resid <- c(.resid, .iniDf$name[.w])
      }
    }
  }
  .est <- .iniDf$name[!is.na(.iniDf$ntheta) & !.iniDf$fix]
  # a covariate coefficient on a theta with an eta is estimated by saem's regression
  .cov <- ui$muRefCovariateDataFrame
  .covOnEta <- .cov$covariateParameter[.cov$theta %in% ui$muRefDataFrame$theta]
  .keep <- c(setdiff(.resid, .predDeps), .llThetas)
  .thetas <- setdiff(intersect(unique(.keep), .est),
                     c(ui$muRefDataFrame$theta, .covOnEta, .skip))
  if (length(.thetas) == 0L) return(.empty)
  .w <- match(.thetas, .iniDf$name)
  .lo <- .iniDf$lower[.w]
  .hi <- .iniDf$upper[.w]
  for (.k in seq_along(.thetas)) {
    .r <- .argRange[[.thetas[.k]]]
    if (is.null(.r)) next
    .lo[.k] <- max(.lo[.k], .r[1])
    .hi[.k] <- min(.hi[.k], .r[2])
  }
  data.frame(theta = .thetas, lower = .lo, upper = .hi)
}

#' Give eta-less thetas a temporary mu-referenced eta on their range's scale
#'
#' Each theta is estimated through an internal `rxBoundedTr.<theta>` and one
#' prepended line on the scale of its range: `exp()` for a positive
#' parameter, `expit()` for a bounded one, additive when unbounded.  A lower
#' bound other than zero takes a helper line, since rxode2 does not
#' mu-reference `a + exp(theta + eta)`.  The transforms are handed to the
#' bounded-transform back-transform (`.postEstimationBoundedTransform()`).
#'
#' @param ui rxode2 ui
#' @param spec data frame from `.saemPseudoEtaThetas()`
#' @param omega initial variance of each temporary eta
#' @return rewritten ui carrying `boundedTransforms`
#' @noRd
.saemAddPseudoEtas <- function(ui, spec, omega = 0.1) {
  .iniDf <- ui$iniDf
  .etaRows <- which(!is.na(.iniDf$neta1))
  .template <- if (length(.etaRows) > 0L) .iniDf[.etaRows[1], , drop = FALSE] else .iniDf[1, , drop = FALSE]
  .maxEta <- if (length(.etaRows) > 0L) max(.iniDf$neta1[.etaRows]) else 0
  .newLines <- character(0)
  .transforms <- vector("list", nrow(spec))
  for (.k in seq_len(nrow(spec))) {
    .t <- spec$theta[.k]
    .lo <- spec$lower[.k]
    .hi <- spec$upper[.k]
    .w <- which(.iniDf$name == .t)
    .est <- .iniDf$est[.w]
    .int <- paste0("rxBoundedTr.", .t)
    .core <- paste0(.int, " + rx.eta.", .t)
    if (is.finite(.lo) && is.finite(.hi)) {
      .eps <- (.hi - .lo) * 1e-6
      .e <- max(.lo + .eps, min(.hi - .eps, .est))
      .type <- "logit"
      .init <- log((.e - .lo) / (.hi - .e))
      .line <- paste0(.t, " <- expit(", .core, ", ", .lo, ", ", .hi, ")")
    } else if (is.finite(.lo)) {
      .type <- "lower_exp"
      .init <- log(max(.est - .lo, 1e-6))
      .line <- if (.lo == 0) {
        paste0(.t, " <- exp(", .core, ")")
      } else {
        c(paste0("rx.l.", .t, " <- exp(", .core, ")"),
          paste0(.t, " <- ", .lo, " + rx.l.", .t))
      }
    } else if (is.finite(.hi)) {
      .type <- "upper_exp"
      .init <- log(max(.hi - .est, 1e-6))
      .line <- paste0(.t, " <- ", .hi, " - exp(", .core, ")")
    } else {
      .type <- "identity"
      .init <- .est
      .line <- paste0(.t, " <- ", .core)
    }
    .newLines <- c(.newLines, .line)
    .transforms[[.k]] <- list(name = .t, internalName = .int, type = .type,
                              lower = .lo, upper = .hi, initTrans = .init,
                              initOrig = .est, pseudoEta = TRUE)
    .iniDf$name[.w] <- .int
    .iniDf$lower[.w] <- -Inf
    .iniDf$upper[.w] <- Inf
    .iniDf$est[.w] <- .init
    .iniDf$err[.w] <- NA_character_
    .iniDf$condition[.w] <- NA_character_
    .maxEta <- .maxEta + 1
    .row <- .template
    .row$ntheta <- NA_integer_
    .row$neta1 <- .row$neta2 <- .maxEta
    .row$name <- paste0("rx.eta.", .t)
    .row$lower <- -Inf
    .row$upper <- Inf
    .row$est <- omega
    .row$fix <- FALSE
    .row$label <- NA_character_
    .row$backTransform <- NA_character_
    .row$condition <- "id"
    .row$err <- NA_character_
    if (any(names(.row) == "prior")) .row$prior <- NA_character_
    .iniDf <- rbind(.iniDf, .row)
  }
  .model <- str2lang(paste0("model({",
                            paste(c(.newLines, vapply(ui$lstExpr, deparse1, character(1))),
                                  collapse = "\n"),
                            "})"))
  .ini <- as.expression(lotri::as.lotri(.iniDf))
  .ini[[1]] <- quote(`ini`)
  .fun <- .getUiFunFromIniAndModel(ui, .ini, .model)
  .newUi <- rxode2::rxUiDecompress(.fun())
  assign("modelName", ui$modelName, envir = .newUi)
  .newUi$boundedTransforms <- .transforms
  .newUi
}

#' Fit a modeled residual error component as its `dnorm()` likelihood in saem
#'
#' The closed-form residual M-step can only estimate a constant residual
#' parameter, so `cp ~ add(a)` with a modeled `a` is rewritten to the
#' equivalent `cp ~ add(a) + dnorm()`.  Any eta-less likelihood theta of a
#' general-likelihood fit then gets a temporary eta (`.saemAddPseudoEtas()`).
#'
#' @param ui rxode2 ui
#' @param est estimation method
#' @param data dataset (unused)
#' @param control control (unused)
#' @return list with the rewritten ui, or `NULL` when nothing changes
#' @noRd
.preProcessSaemModeledResid <- function(ui, est, data, control) {
  if (!identical(est, "saem")) return(NULL)
  .orig <- ui
  .conds <- .saemModeledResidualCond(ui)
  .pred <- ui$predDf
  # dnorm() scores the transformed DV without the lambda-dependent Jacobian
  .lam <- .pred$cond %in% .conds & grepl("boxCox|yeoJohnson", paste(.pred$transform))
  if (any(.lam)) {
    stop("saem cannot fit a modeled residual error with a boxCox()/yeoJohnson() transform ('",
         paste(.pred$cond[.lam], collapse = "', '"), "')", call. = FALSE)
  }
  for (.cond in .conds) {
    .new <- .saemAddDnormToErrLine(ui$lstExpr[[.pred$line[.pred$cond == .cond]]])
    ui <- eval(bquote(rxode2::model(ui, .(.new))))
    warning(sprintf("modeled residual error for '%s'; fit as dnorm() likelihood", .cond),
            call. = FALSE)
  }
  .spec <- .saemPseudoEtaThetas(ui)
  if (length(.conds) == 0L && nrow(.spec) == 0L) return(NULL)
  if (nrow(.spec) > 0L) {
    ui <- .saemAddPseudoEtas(ui, .spec)
    .pre <- "temporary eta for eta-less likelihood theta(s): "
    warning(.pre, .vaeTruncList(.spec$theta, prefix = .pre), call. = FALSE)
  }
  # the reported fit shows the model as written (see .nlmixrEstUpdatesOrigModel)
  if (is.null(nlmixr2global$nlmixr2EstEnv$nlmixrPureInputUi)) {
    nlmixr2global$nlmixr2EstEnv$nlmixrPureInputUi <- rxode2::rxUiDecompress(.orig)
  }
  list(ui = ui)
}

preProcessHooksAdd(".preProcessSaemModeledResid", .preProcessSaemModeledResid)
