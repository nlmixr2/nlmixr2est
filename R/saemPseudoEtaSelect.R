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

#' Record the argument range of a theta that is a likelihood argument itself
#'
#' @param cand candidate list (see `.saemEndpointCandidates()`)
#' @param alias theta the argument is an alias of, or `NA`
#' @param pred1 one `predDf` row
#' @param col the `predDf` column holding the argument
#' @return updated candidate list
#' @noRd
.saemAddArgRange <- function(cand, alias, pred1, col) {
  if (is.na(alias)) return(cand)
  .r <- .saemPredArgRange(pred1, col)
  if (.r$skip) {
    cand$skip <- c(cand$skip, alias)
  } else if (!is.null(.r$range)) {
    cand$argRange[[alias]] <- .r$range
  }
  cand
}

#' Thetas informing the likelihood of one non-`ll()` endpoint
#'
#' @param ui rxode2 ui
#' @param i `predDf` row
#' @param lines deparsed model lines
#' @param lhs `.rxLineLhs(lines)`
#' @return list with `resid` (likelihood thetas), `predDeps` (prediction
#'   dependencies), `skip` and `argRange`
#' @noRd
.saemEndpointCandidates <- function(ui, i, lines, lhs) {
  .pred <- ui$predDf
  .iniDf <- ui$iniDf
  .idx <- .pred$line[i]
  .ret <- list(resid = character(0), skip = character(0), argRange = list(),
               predDeps = .rxMtimeDeps(lines, lhs, .idx, paste(.pred$var[i])))
  for (.c in intersect(c("a", "b", "c", "d", "e", "f"), names(.pred))) {
    .v <- .pred[[.c]][i]
    if (is.na(.v)) next
    .ret$resid <- c(.ret$resid, .rxMtimeDeps(lines, lhs, .idx, .v))
    .ret <- .saemAddArgRange(.ret, .saemAliasTheta(lines, lhs, .idx, .v),
                             .pred[i, , drop = FALSE], .c)
  }
  .err <- which(!is.na(.iniDf$err) & .iniDf$condition == .pred$cond[i])
  .bad <- grepl("boxCox|yeoJohnson|^ar$|^binom$", .iniDf$err[.err])
  .ret$skip <- c(.ret$skip, .iniDf$name[.err[.bad]])
  .ret$resid <- c(.ret$resid, .iniDf$name[.err[!.bad]])
  .ret
}

#' Thetas informing each likelihood endpoint of a model
#'
#' @param ui rxode2 ui
#' @return list with `resid`, `predDeps`, `skip`, `argRange` and `ll`
#' @noRd
.saemLikelihoodCandidates <- function(ui) {
  .pred <- ui$predDf
  .lines <- vapply(ui$lstExpr, deparse1, character(1))
  .lhs <- .rxLineLhs(.lines)
  # an error model `cp ~ add(a)` is not an assignment to cp; treating it as one
  # pulls one endpoint's residual thetas into a later endpoint's prediction
  .lhs[.pred$line] <- NA_character_
  .all <- list(resid = character(0), predDeps = character(0), skip = character(0),
               argRange = list(), ll = character(0))
  for (i in seq_along(.pred$cond)) {
    .dist <- paste(.pred$distribution[i])
    if (.dist == "ordinal") next
    if (.dist == "LL") {
      # ll() has no separate prediction, and a linCmt()/ODE state hides the
      # assignments behind it, so every theta in the model informs it
      .all$ll <- unlist(lapply(ui$lstExpr, all.vars))
      next
    }
    .one <- .saemEndpointCandidates(ui, i, .lines, .lhs)
    for (.n in c("resid", "predDeps", "skip", "argRange")) {
      .all[[.n]] <- c(.all[[.n]], .one[[.n]])
    }
  }
  .all
}

#' Estimated likelihood thetas without an eta in a saem general-likelihood fit
#'
#' saem refines an eta-less theta only through its phi0 step, which leaves a
#' likelihood parameter near its initial value; these thetas get a temporary
#' eta.  A theta that is the likelihood argument itself uses that argument's
#' range; any other theta uses its own `ini()` bounds.  Thetas that feed a
#' prediction (outside `ll()`), count arguments, transform lambdas, ordinal
#' endpoints and covariate coefficients on a theta with an eta are left alone.
#'
#' @param ui rxode2 ui
#' @return data frame with `theta`, `lower` and `upper`
#' @noRd
.saemPseudoEtaThetas <- function(ui) {
  .empty <- data.frame(theta = character(0), lower = numeric(0), upper = numeric(0))
  if (!.saemGeneralLik(ui)) return(.empty)
  .cand <- .saemLikelihoodCandidates(ui)
  .iniDf <- ui$iniDf
  .est <- .iniDf$name[!is.na(.iniDf$ntheta) & !.iniDf$fix]
  # a covariate coefficient on a theta with an eta is estimated by saem's regression
  .cov <- ui$muRefCovariateDataFrame
  .covOnEta <- .cov$covariateParameter[.cov$theta %in% ui$muRefDataFrame$theta]
  .keep <- c(setdiff(.cand$resid, .cand$predDeps), .cand$ll)
  .thetas <- setdiff(intersect(unique(.keep), .est),
                     c(ui$muRefDataFrame$theta, .covOnEta, .cand$skip))
  if (length(.thetas) == 0L) return(.empty)
  .w <- match(.thetas, .iniDf$name)
  .ret <- data.frame(theta = .thetas, lower = .iniDf$lower[.w], upper = .iniDf$upper[.w])
  for (.k in which(.thetas %in% names(.cand$argRange))) {
    .r <- .cand$argRange[[.thetas[.k]]]
    .ret$lower[.k] <- max(.ret$lower[.k], .r[1])
    .ret$upper[.k] <- min(.ret$upper[.k], .r[2])
  }
  .ret
}
