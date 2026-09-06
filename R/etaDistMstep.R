## ODE-free M-step for the parameters of a declared random-effect distribution.
##
## In the (y, eta) augmentation the complete-data likelihood factors as
##
##     log p(y | eta)  +  log p(eta | theta_dist)
##
## and theta_dist appears ONLY in the second term.  Its M-step is therefore a
## pure distribution fit to the sampled etas -- no data term, no ODE solve --
## exactly as saem's residual-error step is a fit to accumulated residuals
## rather than a re-solve.
##
## The `rxEtaDistExpand()` rewrite is what breaks that: writing
## `eta = Q(phi(z))` with `z ~ N(0,1)` moves theta_dist out of the prior and
## into the data likelihood, so it becomes a structural parameter needing an
## ODE solve per objective evaluation and lands in `refinePhi0Lik()`'s
## derivative-free search.  Measured on Bauer's gamma data that search converges
## to a STABLE wrong point (CL 2.49 against a truth of 5.03, unchanged across
## 20/60/200 iterations) and a 12-fold larger evaluation budget does not move it
## -- i.e. the formulation, not the optimizer.
##
## EM lets the augmentation be chosen freely: sample in z-space (good MCMC
## geometry, which is the whole point of Bauer's technique) but take the M-step
## in eta-space.  Both are valid EM algorithms converging to the same MLE.
##
## It is NOT circular.  If z were PRIOR draws then `Q(phi(z))` would be exactly
## `family(theta_old)` and the MLE would return theta_old.  The sampled z are
## POSTERIOR draws -- which is precisely why their spread is measured at ~0.94
## rather than 1.0 -- so the implied etas carry data information.

#' nlmixr2est's own Nelder-Mead, with two guards its interface needs
#'
#' `nmsimplex()` returns a LIST (`$par`, `$value`), and it derives the simplex
#' step as `-0.2 * start` -- so a coordinate whose starting value is 0 gets a
#' step of 0 and can never move.  Offset such coordinates before handing them
#' over.  Used instead of `stats::optim()` so every optimization in this package
#' goes through the same C simplex (neldermead_wrap -> nelder_fn).
#' @noRd
.etaDistNm <- function(start, fn, maxeval = 2000L, reltol = 1e-10) {
  .s <- as.numeric(start)
  .s[!is.finite(.s) | .s == 0] <- 1e-3
  .r <- tryCatch(nmsimplex(.s, fn, control = list(maxeval = maxeval,
                                                  reltol = reltol)),
                 error = function(e) NULL)
  if (is.null(.r) || is.null(.r$par) || anyNA(.r$par) ||
        !all(is.finite(.r$par))) return(NULL)
  list(par = as.numeric(.r$par), value = as.numeric(.r$value))
}

#' Maximum-likelihood M-step for a declared family, from sampled etas
#'
#' Uses the family's OWN density -- `etaDist` records a `d*()` call, so the
#' log-likelihood is that call with `log=TRUE`.  Family-agnostic, and a theta
#' may enter through an arbitrary expression such as `1/exp(lclrv)`.
#'
#' @param etaVals sampled values of the declared random effect
#' @param distCall the `d*()` call recorded on the iniDf
#' @param thetaNames free thetas to estimate
#' @param start their current values
#' @return named vector of updated thetas, or `NULL` if the fit failed
#' @noRd
.etaDistMstep <- function(etaVals, distCall, thetaNames, start) {
  .cl <- if (is.character(distCall)) str2lang(distCall) else distCall
  .dn <- as.character(.cl[[1]])
  if (!exists(.dn, mode = "function")) return(NULL)
  .e <- etaVals[is.finite(etaVals)]
  if (length(.e) < 2L) return(NULL)
  .nll <- function(p) {
    .tv <- stats::setNames(as.list(p), thetaNames)
    .args <- lapply(as.list(.cl)[-1], function(.a) {
      tryCatch(eval(.a, envir = .tv), error = function(e) NA_real_)
    })
    if (any(!vapply(.args, function(.a) is.numeric(.a) && length(.a) == 1L &&
                      is.finite(.a), logical(1)))) return(1e10)
    .ll <- tryCatch(do.call(.dn, c(list(.e), .args, list(log = TRUE))),
                    error = function(e) NULL)
    if (is.null(.ll) || any(!is.finite(.ll))) return(1e10)
    -sum(.ll)
  }
  ## nlmixr2est's own C simplex, matching _saemOpt()/refinePhi0Lik().  This is
  ## the R-side reference implementation used to validate the M-step; the
  ## in-loop version belongs in saem.cpp calling nelder_fn directly, with the
  ## family density coming from Rmath rather than an R callback on the hot path.
  .x <- .etaDistNm(start, .nll)
  if (is.null(.x) || .x$value >= 1e10) return(NULL)
  stats::setNames(.x$par, thetaNames)
}

#' Closed-form M-step for a Gaussian copula's correlation
#'
#' For a copula block the latent pair is bivariate normal with UNIT variances,
#' so the constrained MLE of the correlation is `sum(w1*w2)/n` -- a closed form,
#' not a search.  This is the same quantity the plan specifies as a simulation
#' CHECK (`cor(qnorm(pgamma(cl)), qnorm(pgamma(v1)))`); it works equally well as
#' an estimator.
#' @noRd
.etaDistCorMstep <- function(w1, w2) {
  .ok <- is.finite(w1) & is.finite(w2)
  if (sum(.ok) < 2L) return(NULL)
  .r <- sum(w1[.ok] * w2[.ok]) / sum(.ok)
  if (!is.finite(.r)) return(NULL)
  max(min(.r, 0.999), -0.999)
}

#' Family code for the C++ distribution M-step
#'
#' The code IS the row number in `lotri::lotriEtaDists()`, so the C++ dispatch
#' (`RXETADIST_*`, src/saem.cpp) and the catalog cannot drift apart -- adding a
#' family to the catalog shifts nothing already assigned.
#'
#' Returns `0L` for anything the C++ dispatch does not implement, which selects
#' the general R fallback (`.etaDistMstep()`).  That fallback evaluates the
#' declaration's own `d*()` call, so a family is never LESS supported than it
#' was; it just pays an R round trip per objective evaluation.
#' @noRd
.etaDistFamilyCode <- function(distCall) {
  .cl <- if (is.character(distCall)) str2lang(distCall) else distCall
  .nm <- as.character(.cl[[1]])
  .tab <- try(lotri::lotriEtaDists(), silent = TRUE)
  if (inherits(.tab, "try-error") || is.null(.tab$name)) return(0L)
  .w <- which(.tab$name == .nm)
  if (length(.w) != 1L) return(0L)
  as.integer(.w)
}

#' Support of a declared family ("real", "positive", "nonneg", "unit")
#'
#' Used to pick the surrogate in [etaDistInit()]: a positive-support family
#' gets a log-normal surrogate, a real-support one an ordinary normal.
#' @noRd
.etaDistSupport <- function(distCall) {
  .cl <- if (is.character(distCall)) str2lang(distCall) else distCall
  .nm <- as.character(.cl[[1]])
  .tab <- try(lotri::lotriEtaDists(), silent = TRUE)
  if (inherits(.tab, "try-error") || is.null(.tab$support)) return(NA_character_)
  .w <- which(.tab$name == .nm)
  if (length(.w) != 1L) return(NA_character_)
  .tab$support[.w]
}

#' Solve a declared family's thetas so its ARGUMENTS take given values
#'
#' The C++ M-step estimates the family's NATIVE parameters (shape, rate, ...);
#' this turns those back into the user's thetas, which may enter through
#' arbitrary expressions (`shape = 1/exp(lclrv)`).  Called ONCE per iteration,
#' not per objective evaluation -- the expensive part (the likelihood over
#' every sampled eta) stays in C++.
#'
#' Exact when the map is invertible; when it is not, this returns the
#' least-squares closest thetas, which is the best that can be done without
#' constraining what `dist()` accepts.
#' @noRd
.etaDistArgsToThetas <- function(distCall, thetaNames, start, targetArgs) {
  .cl <- if (is.character(distCall)) str2lang(distCall) else distCall
  .ex <- as.list(.cl)[-1]
  if (length(.ex) != length(targetArgs)) return(NULL)
  .obj <- function(p) {
    .tv <- stats::setNames(as.list(p), thetaNames)
    .a <- vapply(.ex, function(.e) {
      .v <- tryCatch(eval(.e, envir = .tv), error = function(e) NA_real_)
      if (!is.numeric(.v) || length(.v) != 1L) NA_real_ else as.numeric(.v)
    }, numeric(1))
    if (anyNA(.a) || any(!is.finite(.a))) return(1e10)
    ## relative, so arguments on very different scales weigh comparably
    sum((log(pmax(abs(.a), 1e-300)) - log(pmax(abs(targetArgs), 1e-300)))^2) +
      sum((sign(.a) != sign(targetArgs)) * 1e3)
  }
  .x <- .etaDistNm(start, .obj)
  if (is.null(.x) || .x$value > 1e-6) return(NULL)
  stats::setNames(.x$par, thetaNames)
}

#' Metadata the C++ distribution M-step needs
#'
#' Returns the phi column of each declared eta's own latent normal, its family
#' code, which declared eta it is copula-correlated with, its current NATIVE
#' parameters and correlation -- plus the bookkeeping R needs to map the
#' updated parameters back onto thetas.
#' @noRd
.etaDistMstepInfo <- function(ui, etaTrans, etaNames) {
  .ui <- rxode2::rxUiDecompress(ui)
  .d <- rxode2::rxUiEtaDists(.ui)
  if (nrow(.d) == 0L) return(NULL)
  .ini <- .ui$iniDf
  .thNames <- .ini$name[!is.na(.ini$ntheta)]
  .thVals <- stats::setNames(as.list(.ini$est[!is.na(.ini$ntheta)]), .thNames)
  .n <- nrow(.d)
  .lat <- integer(.n); .fam <- integer(.n); .cw <- rep(-1L, .n)
  .maxA <- 0L; .args <- vector("list", .n); .tn <- vector("list", .n)
  for (.i in seq_len(.n)) {
    .e <- .d$name[.i]
    ## the expansion renames the declared eta's latent `rxz.<eta>`
    .lz <- paste0("rxz.", .e)
    .w <- which(etaNames == .lz)
    .lat[.i] <- if (length(.w) == 1L) as.integer(etaTrans[.w]) - 1L else -1L
    .fam[.i] <- .etaDistFamilyCode(.d$etaDist[.i])
    .cl <- str2lang(.d$etaDist[.i])
    .a <- vapply(as.list(.cl)[-1], function(.x) {
      .v <- tryCatch(eval(.x, envir = .thVals), error = function(e) NA_real_)
      if (!is.numeric(.v) || length(.v) != 1L) NA_real_ else as.numeric(.v)
    }, numeric(1))
    .args[[.i]] <- .a
    .maxA <- max(.maxA, length(.a))
    .tn[[.i]] <- intersect(all.vars(.cl), .thNames)
  }
  if (any(.lat < 0) || any(.fam <= 0)) return(NULL)   # -> R fallback
  ## Copula structure.  In the DECLARED model a copula block is an omega block
  ## with a unit diagonal whose off-diagonal IS the correlation; the expansion
  ## turns that into independent latents plus an in-model Cholesky, so the
  ## correlation has to be carried here for the driver to rebuild the
  ## correlated latent.  Only a PAIR is handled -- etaDistMstep() reconstructs
  ## a single partner, so a larger block would be silently wrong; return NULL
  ## and let the R fallback take it.
  .rho <- rep(0, .n)
  .off <- .ini[!is.na(.ini$neta1) & !is.na(.ini$neta2) &
                 .ini$neta1 != .ini$neta2, , drop = FALSE]
  .netaOf <- function(.nm) {
    .w <- which(.ini$name == .nm & .ini$neta1 == .ini$neta2)
    if (length(.w) == 1L) as.integer(.ini$neta1[.w]) else NA_integer_
  }
  .id <- vapply(.d$name, .netaOf, integer(1))
  if (nrow(.off) > 0L) {
    for (.r in seq_len(nrow(.off))) {
      .a1 <- which(.id == .off$neta1[.r]); .a2 <- which(.id == .off$neta2[.r])
      if (length(.a1) != 1L || length(.a2) != 1L) next   # not a declared pair
      .hi <- max(.a1, .a2); .lo <- min(.a1, .a2)
      if (.cw[.hi] >= 0L) return(NULL)                   # >2 declared partners
      .cw[.hi] <- .lo - 1L                               # 0-based
      .rho[.hi] <- .off$est[.r]
    }
  }
  .am <- matrix(0, nrow = .n, ncol = max(1L, .maxA))
  for (.i in seq_len(.n)) .am[.i, seq_along(.args[[.i]])] <- .args[[.i]]
  list(latent = .lat, fam = .fam, corWith = .cw,
       args = .am, rho = .rho,
       dist = .d$etaDist, thetas = .tn, etas = .d$name)
}
