## Two-stage initialization for declared non-normal random effects.
##
## The estimation problem these models present is almost entirely a starting
## value problem for the DISTRIBUTION parameters.  Measured on Bauer's gamma
## dataset: started at NONMEM's own estimates the latent normals come out
## correctly dispersed (sd 0.94 against NONMEM's 0.94-0.99) and the copula
## correlation lands at 0.459 against a truth of 0.438, while a cold start
## overshoots the correlation to 0.78-0.85 and leaves CL near its ini value.
##
## So: solve a cheap surrogate first and hand its answer over.  The surrogate
## replaces each declared family with a LOG-NORMAL, which `rxEtaDistExpand()`
## collapses to `exp(mu + sd*z)` -- no inverse CDF, no `phiU()` round trip --
## so it is an ordinary, well-conditioned saem fit through machinery that is
## already well tested.  The copula, the model rewriting and everything else
## stay exactly as they are.
##
## Deliberately NOT posthoc/MAP: empirical Bayes estimates are shrunk toward
## the very distribution being estimated, so fitting a family to them would
## inherit that shrinkage and understate the relative variance.  A surrogate
## POPULATION fit has no such circularity.

#' Gauss-Hermite nodes/weights for standard-normal expectations
#'
#' `E[f(Z)] = sum(w * f(z))` for `Z ~ N(0,1)`.  Same scaling the AGQ code uses
#' (`.agq()`, R/agq.R): fastGHQuad's nodes are for `exp(-x^2)`, so they scale by
#' `sqrt(2)` and the weights by `1/sqrt(pi)`.
#' @noRd
.etaDistGh <- function(n = 30L) {
  rxode2::rxReq("fastGHQuad")
  .gh <- fastGHQuad::gaussHermiteData(n)
  list(z = .gh$x * sqrt(2), w = .gh$w / sqrt(pi))
}

#' Mean and variance of a declared family, by quadrature over the latent normal
#'
#' `etaDist` is the `d*()` call recorded on the iniDf; its arguments are
#' expressions over the model's thetas.  The corresponding quantile function is
#' the base-R `q*()` sibling, so `E[eta]` and `Var[eta]` follow from
#' `eta = Q(pnorm(z))` integrated against the standard normal.
#' @noRd
.etaDistMoments <- function(distCall, thetaVals, gh = .etaDistGh()) {
  .cl <- if (is.character(distCall)) str2lang(distCall) else distCall
  .qn <- sub("^d", "q", as.character(.cl[[1]]))
  if (!exists(.qn, mode = "function")) return(NULL)
  .args <- lapply(as.list(.cl)[-1], function(.a) {
    tryCatch(eval(.a, envir = as.list(thetaVals)), error = function(e) NA_real_)
  })
  if (any(!vapply(.args, function(.a) is.numeric(.a) && length(.a) == 1L &&
                    is.finite(.a), logical(1)))) return(NULL)
  ## The outermost Gauss-Hermite node saturates pnorm() to exactly 1 in double
  ## precision, and an inverse CDF at 1 is Inf -- the same boundary Bauer
  ## guards with DEL and rxode2 with phiU()'s clamp.  Clamp identically
  ## (RX_PHIU_EPS); the affected node carries weight ~1e-22, so it changes
  ## nothing except the overflow.
  .eps <- 1e-15
  .u <- pmin(pmax(stats::pnorm(gh$z), .eps), 1 - .eps)
  .q <- tryCatch(do.call(.qn, c(list(.u), .args)), error = function(e) NULL)
  if (is.null(.q) || !all(is.finite(.q))) return(NULL)
  .m <- sum(gh$w * .q)
  .v <- sum(gh$w * (.q - .m)^2)
  ## LOG-scale moments too, and they are the ones the surrogate should use.
  ##
  ## The surrogate IS a log-normal, so meanlog/varlog ARE its parameters -- on
  ## that side the conversion is an identity rather than an approximation.  The
  ## arithmetic pair is not interchangeable with it: a log-normal and a gamma
  ## sharing a log-variance have wildly different relative variances once the
  ## family is heavy-tailed.  Measured on Bauer's g4 (gamma shape 0.5, relvar
  ## 2.0), matching the arithmetic pair returns relvar 68 where matching the
  ## log pair returns 1.83.
  ##
  ## NA for a family whose support reaches 0 or below, where log is undefined;
  ## callers fall back to the arithmetic pair.
  .lm <- NA_real_; .lv <- NA_real_
  if (all(.q > 0)) {
    .lq <- log(.q)
    .lm <- sum(gh$w * .lq)
    .lv <- sum(gh$w * (.lq - .lm)^2)
  }
  c(mean = .m, var = .v, meanlog = .lm, varlog = .lv)
}

#' Solve a declared family's thetas so its mean/variance match a target
#'
#' Generic across families: the family's own quantile function supplies the
#' moments (see .etaDistMoments), so no per-family algebra is needed and a
#' theta may enter through an arbitrary expression such as `1/exp(lclrv)`.
#' @noRd
.etaDistSolveThetas <- function(distCall, thetaNames, start, target,
                                gh = .etaDistGh()) {
  .useLog <- all(c("meanlog", "varlog") %in% names(target)) &&
    all(is.finite(target[c("meanlog", "varlog")]))
  .obj <- function(p) {
    .tv <- stats::setNames(as.list(p), thetaNames)
    .mv <- .etaDistMoments(distCall, .tv, gh)
    if (is.null(.mv)) return(1e10)
    if (.useLog) {
      if (!all(is.finite(.mv[c("meanlog", "varlog")]))) return(1e10)
      ## meanlog is a location on the log scale and can be zero or negative, so
      ## it is matched by DIFFERENCE; varlog is positive and matched by ratio,
      ## which keeps both terms dimensionless the way the arithmetic pair was.
      (.mv[["meanlog"]] - target[["meanlog"]])^2 +
        (log(.mv[["varlog"]] / target[["varlog"]]))^2
    } else {
      ## relative, so a mean of 5 and a variance of 2 weigh comparably
      (log(.mv[["mean"]] / target[["mean"]]))^2 +
        (log(.mv[["var"]] / target[["var"]]))^2
    }
  }
  ## nlmixr2est's own Nelder-Mead (nmsimplex -> neldermead_wrap -> nelder_fn),
  ## not stats::optim: every other optimization in this package goes through
  ## the C simplex and this should not be the one exception.
  .x <- .etaDistNm(start, .obj)
  if (is.null(.x) || .x$value > 1e-3) return(NULL)
  stats::setNames(.x$par, thetaNames)
}

#' Build the ordinary-random-effect surrogate of a declared model
#'
#' Each declared eta becomes a PLAIN random effect: the `dist()` declaration is
#' dropped, the parameter is rewritten `p <- exp(mu + eta)`, and the latent
#' block's fixed unit diagonal is freed.  `mu` is then mu-referenced (theta+eta,
#' so saem estimates it with the closed-form GLS) and the spread comes from the
#' ordinary omega M-step -- both well-conditioned closed forms.
#'
#' Substituting `dist(eta) ~ dlnorm(...)` instead does NOT work, and it is worth
#' recording why: a log-normal declaration collapses the inverse-CDF transform,
#' but its parameters remain non-mu thetas inside a declaration, so they route
#' to phi0 and `refinePhi0Lik()` -- the same derivative-free search that
#' mis-estimates the declared model.  Measured on Bauer's gamma data, that
#' surrogate returned a log-normal sd of 0.067 where ~0.29 was wanted, i.e. it
#' inherited the very problem it was meant to sidestep.
#'
#' Only the "the parameter IS the variate" shape is handled (`p <- eta`, which
#' is what a positive-support declaration such as `dgamma()` produces).  A
#' declaration used as `p <- exp(theta + eta)` already has a normal surrogate --
#' just drop the declaration -- and is returned unchanged.
#' @noRd
.etaDistSurrogate <- function(ui, gh = .etaDistGh()) {
  .ui <- rxode2::rxUiDecompress(ui)
  .d <- rxode2::rxUiEtaDists(.ui)
  if (nrow(.d) == 0L) return(NULL)
  .ini <- .ui$iniDf
  .thNames <- .ini$name[!is.na(.ini$ntheta)]
  .thVals <- stats::setNames(as.list(.ini$est[!is.na(.ini$ntheta)]), .thNames)
  .iniTxt <- deparse(.ui$iniFun)
  .modTxt <- deparse(.ui$modelFun)
  .seed <- list()
  for (.i in seq_len(nrow(.d))) {
    .e <- .d$name[.i]
    .mv <- .etaDistMoments(.d$etaDist[.i], .thVals, gh)
    if (is.null(.mv) || !is.finite(.mv[["mean"]]) ||
          !is.finite(.mv[["var"]]) || .mv[["mean"]] <= 0) return(NULL)
    ## Seed the log-normal from the declared family's LOG-scale moments, which
    ## are its parameters outright.  Going through the arithmetic pair instead
    ## understates the log-spread badly for a heavy tail: on g4's ini() it gives
    ## varlog 1.10 where the declared gamma's own varlog is 4.97.
    if (all(is.finite(.mv[c("meanlog", "varlog")])) && .mv[["varlog"]] > 0) {
      .mu <- .mv[["meanlog"]]; .sd <- sqrt(.mv[["varlog"]])
    } else {
      .cv2 <- .mv[["var"]] / (.mv[["mean"]]^2)
      .sd <- sqrt(log1p(.cv2)); .mu <- log(.mv[["mean"]]) - 0.5 * .sd^2
    }
    ## the model line whose entire right-hand side is this eta
    .pat <- paste0("^(\\s*[A-Za-z._][A-Za-z0-9._]*\\s*<-\\s*)",
                   gsub("[.]", "[.]", .e), "\\s*$")
    .w <- grep(.pat, .modTxt)
    if (length(.w) != 1L) return(NULL)
    .mn <- paste0("rxWs.mu.", .e)
    .modTxt[.w] <- sub(.pat, paste0("\\1exp(", .mn, " + ", .e, ")"), .modTxt[.w])
    .seed[[.e]] <- c(mu = .mn, sd = .sd)
    .iniTxt <- c(.iniTxt[-length(.iniTxt)],
                 paste0("  ", .mn, " <- ", format(.mu, digits = 10)), "})")
    ## drop the declaration
    .iniTxt <- .iniTxt[!grepl(paste0("^\\s*dist\\(\\s*", gsub("[.]", "[.]", .e),
                                     "\\s*\\)\\s*~"), .iniTxt)]
  }
  ## Free the latent block.  Its unit diagonal exists only so the copula is a
  ## correlation matrix; as an ordinary random effect it is estimated, and its
  ## starting value should be the spread the declared family actually implies
  ## (from .seed) rather than 1 -- otherwise the block starts at a correlation
  ## near 1 and the surrogate begins somewhere the declared model never was.
  for (.e in .d$name) {
    .sd <- .seed[[.e]][["sd"]]
    .p <- paste0("^(\\s*", gsub("[.]", "[.]", .e), "\\s*(\\+[^~]*)?~\\s*)")
    .k <- grep(.p, .iniTxt)
    if (length(.k) != 1L) next
    ## Under the unit diagonal the off-diagonals ARE correlations, so freeing
    ## the block means rescaling: diagonal -> sd^2, off-diagonal -> rho*sd_i*sd_j.
    ## Left as correlations they exceed the new variances and the block is not
    ## positive definite.
    .sdj <- as.numeric(.sd)
    .num <- regmatches(.iniTxt[.k], regexpr("~.*$", .iniTxt[.k]))
    .vals <- suppressWarnings(as.numeric(strsplit(gsub("[~c()\\s]", "", .num),
                                                  ",")[[1]]))
    if (anyNA(.vals)) next
    ## row j of the block: j-1 correlations against the earlier declared etas,
    ## then this eta's own (unit) variance
    .prior <- .d$name[seq_len(length(.vals) - 1L)]
    .sdi <- vapply(.prior, function(.q) as.numeric(.seed[[.q]][["sd"]]), numeric(1))
    .vals[seq_len(length(.vals) - 1L)] <- .vals[seq_len(length(.vals) - 1L)] * .sdi * .sdj
    .vals[length(.vals)] <- .sdj^2
    .iniTxt[.k] <- sub("~.*$",
                       if (length(.vals) == 1L) paste0("~ ", format(.vals, digits = 10))
                       else paste0("~ c(", paste(format(.vals, digits = 10),
                                                 collapse = ", "), ")"),
                       .iniTxt[.k])
  }
  ## thetas orphaned by dropping the declarations
  .cand <- intersect(unique(unlist(lapply(.d$etaDist, function(.s) all.vars(str2lang(.s))),
                                   use.names = FALSE)), .thNames)
  for (.c in .cand) {
    .isIni <- grepl(paste0("^\\s*", gsub("[.]", "[.]", .c), "\\s*(<-|~)"), .iniTxt)
    if (!any(grepl(paste0("\\b", .c, "\\b"), c(.iniTxt[!.isIni], .modTxt)))) {
      .iniTxt <- .iniTxt[!.isIni]
    }
  }
  .fun <- try(eval(parse(text = paste0("function() {\n",
                                       paste(.iniTxt, collapse = "\n"), "\n",
                                       paste(.modTxt, collapse = "\n"), "\n}"))),
              silent = TRUE)
  if (inherits(.fun, "try-error")) return(NULL)
  list(fun = .fun, seed = .seed)
}

#' Two-stage starting values for a declared non-normal random effect
#'
#' Fits a log-normal surrogate of `object` and moment-matches its answer back
#' onto the declared families, returning a model with better `ini()` values.
#' Intended as a warm start for `est="saem"`/`"focei"` on `dist()` models,
#' whose cold-start difficulty is largely a starting-value problem.
#'
#' @param object model with at least one `dist()` declaration
#' @param data data to fit the surrogate against
#' @param control control for the surrogate fit; the defaults are deliberately
#'   short, since only starting values are wanted
#' @param ... passed to [nlmixr2()]
#' @return `object` with updated `ini()` estimates
#' @export
etaDistInit <- function(object, data, control = saemControl(nBurn = 100, nEm = 100,
                                                            print = 0, covMethod = ""),
                        ...) {
  .ui <- rxode2::rxUiDecompress(rxode2::assertRxUi(object))
  .d <- rxode2::rxUiEtaDists(.ui)
  if (nrow(.d) == 0L) {
    stop("'etaDistInit()' needs a model with at least one 'dist()' declaration",
         call. = FALSE)
  }
  .gh <- .etaDistGh()
  .s <- .etaDistSurrogate(.ui, .gh)
  if (is.null(.s)) {
    warning("could not build a log-normal surrogate; starting values unchanged",
            call. = FALSE)
    return(object)
  }
  .fit <- try(suppressWarnings(nlmixr2(.s$fun, data, est = "saem",
                                       control = control, ...)), silent = TRUE)
  if (inherits(.fit, "try-error")) {
    warning("the surrogate fit failed; starting values unchanged\n  ",
            conditionMessage(attr(.fit, "condition")), call. = FALSE)
    return(object)
  }
  .est <- .fit$parFixedDf[, "Estimate"]
  names(.est) <- rownames(.fit$parFixedDf)
  .ini <- .ui$iniDf
  .thNames <- .ini$name[!is.na(.ini$ntheta)]
  for (.i in seq_len(nrow(.d))) {
    .e <- .d$name[.i]
    .mn <- .s$seed[[.e]][["mu"]]
    if (!(.mn %in% names(.est))) next
    .mu <- .est[[.mn]]
    ## The surrogate's spread is an ORDINARY omega now, not a theta.
    .om <- .fit$omega
    if (is.null(.om) || !(.e %in% rownames(.om))) next
    .w <- .om[.e, .e]
    if (!is.finite(.w) || .w <= 0) next
    ## The surrogate's fitted meanlog/varlog ARE (.mu, .w) -- pass them across
    ## unchanged and let the declared family match them on the log scale.  The
    ## arithmetic pair is carried too, as the fallback for a family whose
    ## support reaches zero, but it must not be the primary target: exp(.w) - 1
    ## turns the surrogate's varlog 4.24 into relvar 68 where the declared
    ## gamma's answer is 1.83 against a truth of 2.0.
    .mean <- exp(.mu + 0.5 * .w)
    .var <- (exp(.w) - 1) * exp(2 * .mu + .w)
    if (!is.finite(.mean) || !is.finite(.var) || .var <= 0) next
    .tn <- intersect(all.vars(str2lang(.d$etaDist[.i])), .thNames)
    if (length(.tn) == 0L) next
    .start <- stats::setNames(.ini$est[match(.tn, .ini$name)], .tn)
    .sol <- .etaDistSolveThetas(.d$etaDist[.i], .tn, .start,
                                c(mean = .mean, var = .var,
                                  meanlog = .mu, varlog = .w), .gh)
    if (is.null(.sol)) {
      warning("could not match '", .e, "' to the surrogate's moments; ",
              "its starting values are unchanged", call. = FALSE)
      next
    }
    ## Keep the answer only if it is actually CLOSER to the surrogate than what
    ## the model already had.  A warm start has one job, and there is no reason
    ## to accept one that fails at it -- silently, as this did: the surrogate
    ## fit succeeds, the solver reports convergence because it matched the
    ## moments it was handed, and nothing downstream is placed to notice the
    ## values are absurd.  Before this check the g4 arm was seeded with a
    ## relative variance of 68 against a truth of 2, and the fit it warm-started
    ## ran away to CL 2312 (MARE 9132%) where the COLD start reached 25.9%.
    ##
    ## Threshold-free on purpose: the comparison is against the model's own
    ## starting values on the surrogate's own criterion, so it needs no notion
    ## of what counts as an implausible parameter for an arbitrary family.
    .fitTo <- function(.v) {
      .mm <- .etaDistMoments(.d$etaDist[.i], as.list(.v), .gh)
      if (is.null(.mm) || !all(is.finite(.mm[c("meanlog", "varlog")])) ||
            .mm[["varlog"]] <= 0) return(Inf)
      (.mm[["meanlog"]] - .mu)^2 + (log(.mm[["varlog"]] / .w))^2
    }
    if (all(is.finite(c(.mu, .w))) && .w > 0 &&
          !(.fitTo(.sol) < .fitTo(.start))) {
      warning("the surrogate's starting values for '", .e, "' were no better ",
              "than the model's own; they are unchanged", call. = FALSE)
      next
    }
    for (.t in .tn) .ini$est[.ini$name == .t] <- .sol[[.t]]
  }
  ## The copula correlation of a Gaussian copula IS the correlation of the
  ## latent normals, and for a log-normal surrogate the log-scale correlation is
  ## exactly that -- so it reads straight off the surrogate's fitted omega, with
  ## no moment matching.  The declared block stores it on the unit-diagonal
  ## (correlation) scale.
  .om <- .fit$omega
  .en <- .d$name
  if (!is.null(.om) && length(.en) > 1L && all(.en %in% rownames(.om))) {
    for (.a in seq_len(length(.en) - 1L)) {
      for (.b in seq(.a + 1L, length(.en))) {
        .r <- .om[.en[.a], .en[.b]] /
          sqrt(.om[.en[.a], .en[.a]] * .om[.en[.b], .en[.b]])
        if (!is.finite(.r)) next
        .r <- max(min(.r, 0.999), -0.999)
        .w <- which(.ini$name == .en[.b] &
                      .ini$neta1 != .ini$neta2)
        .w2 <- which(!is.na(.ini$neta1) & !is.na(.ini$neta2) &
                       .ini$neta1 != .ini$neta2 &
                       ((.ini$name == paste0("(", .en[.a], ",", .en[.b], ")")) |
                          (.ini$name == paste0("(", .en[.b], ",", .en[.a], ")"))))
        if (length(.w2) == 1L) .ini$est[.w2] <- .r
      }
    }
  }
  assign("iniDf", .ini, envir = .ui)
  .ui
}
