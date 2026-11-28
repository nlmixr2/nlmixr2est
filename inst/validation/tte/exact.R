# Exact marginal maximum likelihood for the tte validation scenarios.
#
# An implementation independent of rxode2 and nlmixr2est: each scenario's
# cumulative hazard in closed form, integrated over its one random effect by
# Gauss-Hermite quadrature.  `pkJoint` (three random effects) has none.

# Gauss-Hermite nodes and weights for a standard normal (Golub-Welsch)
.tteGaussHermite <- function(n = 40L) {
  .J <- matrix(0, n, n)
  for (.i in seq_len(n - 1L)) {
    .J[.i, .i + 1L] <- .J[.i + 1L, .i] <- sqrt(.i)
  }
  .e <- eigen(.J, symmetric = TRUE)
  list(x = .e$values, w = .e$vectors[1, ]^2)
}

# Log-likelihood of each record given the cumulative hazard H(t) and hazard
# h(t) at the record times of one endpoint.  The first record starts the
# observation and scores nothing.  `H` and `h` are matrices: records x nodes.
.tteRecordLl <- function(dv, H, h, interval = FALSE) {
  .k <- nrow(H)
  if (.k < 2L) {
    return(rep(0, ncol(H)))
  }
  .dH <- H[-1, , drop = FALSE] - H[-.k, , drop = FALSE]
  .dv <- dv[-1]
  .ll <- -.dH
  if (any(.dv == 1)) {
    .e <- which(.dv == 1)
    .ll[.e, ] <- if (interval) {
      log(-expm1(-.dH[.e, , drop = FALSE]))
    } else {
      log(h[-1, , drop = FALSE][.e, , drop = FALSE]) - .dH[.e, , drop = FALSE]
    }
  }
  colSums(.ll)
}

# Weibull cumulative hazard and hazard at times t for scales sc (one per node)
.tteWeibullHh <- function(t, sc, sh) {
  .u <- outer(t, sc, "/")
  list(H = .u^sh, h = sweep(sh * .u^(sh - 1), 2, sc, "/"))
}

# per scenario: parameter names (theta..., then the eta variance) and the
# conditional log-likelihood of one subject at each node value of the eta
.tteExactModels <- list(
  weibull = list(
    par = c("lsc", "lsh", "eta.sc"),
    ll = function(p, s, eta) {
      .w <- .tteWeibullHh(s$TIME, exp(p[1] + eta), exp(p[2]))
      .tteRecordLl(s$DV, .w$H, .w$h)
    }
  ),
  lateEntry = list(
    par = c("lsc", "lsh", "eta.sc"),
    ll = function(p, s, eta) {
      .w <- .tteWeibullHh(s$TIME, exp(p[1] + eta), exp(p[2]))
      .tteRecordLl(s$DV, .w$H, .w$h)
    }
  ),
  interval = list(
    par = c("lsc", "lsh", "eta.sc"),
    ll = function(p, s, eta) {
      .w <- .tteWeibullHh(s$TIME, exp(p[1] + eta), exp(p[2]))
      .tteRecordLl(s$DV, .w$H, .w$h, interval = TRUE)
    }
  ),
  rtte = list(
    par = c("lh", "eta.h"),
    ll = function(p, s, eta) {
      .lam <- exp(p[1] + eta)
      .tteRecordLl(s$DV, outer(s$TIME, .lam), matrix(.lam, nrow(s), length(eta), byrow = TRUE))
    }
  ),
  gapTime = list(
    par = c("lte", "lsh", "eta.te"),
    ll = function(p, s, eta) {
      .te <- exp(p[1] + eta)
      .sh <- exp(p[2])
      # each interval (t[k-1], t[k]] runs on the clock started at TLAST[k]
      .k <- nrow(s)
      if (.k < 2L) {
        return(rep(0, length(eta)))
      }
      .lo <- outer(s$TIME[-.k] - s$TLAST[-1], .te, "/")
      .hi <- outer(s$TIME[-1] - s$TLAST[-1], .te, "/")
      .dH <- .hi^.sh - .lo^.sh
      .ll <- -.dH
      .e <- which(s$DV[-1] == 1)
      .ll[.e, ] <- .ll[.e, ] + log(sweep(.sh * .hi[.e, , drop = FALSE]^(.sh - 1), 2, .te, "/"))
      colSums(.ll)
    }
  ),
  pkFixed = list(
    par = c("llam", "lic50", "eta.lam"),
    ll = function(p, s, eta) {
      .ic50 <- exp(p[2])
      .k <- s$CLI[1] / s$VI[1]
      .c0 <- 100 / s$VI[1]
      .c <- .c0 * exp(-.k * s$TIME)
      # the integral of 1 - C/(IC50 + C) for a mono-exponential C
      .int <- s$TIME - log((.ic50 + .c0) / (.ic50 + .c)) / .k
      .lam <- exp(p[1] + eta)
      .tteRecordLl(s$DV, outer(.int, .lam), outer(1 - .c / (.ic50 + .c), .lam))
    }
  ),
  tvCov = list(
    par = c("lh", "bx", "eta.h"),
    ll = function(p, s, eta) {
      .tsw <- s$TSW[1]
      .int <- pmin(s$TIME, .tsw) + exp(p[2]) * pmax(s$TIME - .tsw, 0)
      .lam <- exp(p[1] + eta)
      # the hazard just before each record
      .rel <- exp(p[2] * (s$TIME > .tsw))
      .tteRecordLl(s$DV, outer(.int, .lam), outer(.rel, .lam))
    }
  ),
  dropout = list(
    par = c("lh1", "lh2", "eta.h1"),
    ll = function(p, s, eta) {
      .s1 <- s[s$DVID == 1, ]
      .s2 <- s[s$DVID == 2, ]
      .lam <- exp(p[1] + eta)
      .l1 <- .tteRecordLl(.s1$DV, outer(.s1$TIME, .lam), matrix(.lam, nrow(.s1), length(eta), byrow = TRUE))
      .l2 <- .tteRecordLl(.s2$DV, matrix(.s2$TIME * exp(p[2])), matrix(exp(p[2]), nrow(.s2)))
      .l1 + .l2
    }
  )
)

#' Exact marginal maximum likelihood fit of a validation scenario
#'
#' @param name scenario name
#' @param dir the validation directory (holding `data/`)
#' @param fixed named parameter values held fixed (e.g. `c(bx = 0)`)
#' @param start starting values (named; the truth when NULL)
#' @return a list with `est` (thetas and the eta variance), `se` and `m2ll`
#'   (-2 log-likelihood); NULL when the scenario has no exact likelihood
#' @noRd
tteValidationExact <- function(name, dir, fixed = NULL, start = NULL) {
  .d <- utils::read.csv(file.path(dir, "data", paste0(name, ".csv")))
  .d <- .d[.d$EVID == 0, ]
  if (name == "exponential") {
    # events / time at risk
    .s <- split(.d, .d$ID)
    .ev <- sum(.d$DV)
    .risk <- sum(vapply(.s, function(s) max(s$TIME) - min(s$TIME), numeric(1)))
    .lh <- log(.ev / .risk)
    return(list(
      est = c(lh = .lh),
      se = c(lh = 1 / sqrt(.ev)),
      m2ll = -2 * (.ev * .lh - exp(.lh) * .risk)
    ))
  }
  .m <- .tteExactModels[[name]]
  if (is.null(.m)) {
    return(NULL)
  }
  .gh <- .tteGaussHermite()
  .subj <- split(.d, factor(.d$ID, levels = unique(.d$ID)))
  if (name != "dropout") {
    .subj <- lapply(.subj, function(s) s[s$DVID == 1, ])
  }
  .np <- length(.m$par)
  # optimize thetas and log(omega)
  .free <- setdiff(.m$par, names(fixed))
  .full <- function(q) {
    .p <- stats::setNames(numeric(.np), .m$par)
    .p[.free] <- q
    .p[names(fixed)] <- fixed
    # the eta variance is optimized on the log scale
    if (.m$par[.np] %in% names(fixed)) {
      .p[.np] <- log(.p[.np])
    }
    .p
  }
  .nll <- function(q) {
    .p <- .full(q)
    .sd <- exp(.p[.np] / 2)
    .eta <- .gh$x * .sd
    -sum(vapply(
      .subj,
      function(s) {
        .l <- .m$ll(.p, s, .eta)
        .mx <- max(.l)
        .mx + log(sum(.gh$w * exp(.l - .mx)))
      },
      numeric(1)
    ))
  }
  if (is.null(start)) {
    .sc <- tteValidationScenarios(name)[[1]]
    start <- rxode2::rxode2(.sc$model)$iniDf
    start <- stats::setNames(start$est, start$name)
  }
  .q0 <- start[.free]
  .q0[.free == .m$par[.np]] <- log(.q0[.free == .m$par[.np]])
  if (length(.free) == 0L) {
    # every parameter fixed: the -2LL at those values
    return(list(est = fixed[.m$par], se = stats::setNames(rep(NA_real_, .np), .m$par), m2ll = 2 * .nll(numeric(0))))
  }
  .o <- stats::optim(.q0, .nll, method = "BFGS", control = list(reltol = 1e-12, maxit = 1000))
  .H <- stats::optimHess(.o$par, .nll)
  .se <- sqrt(diag(solve(.H)))
  .est <- .full(.o$par)
  .est[.np] <- exp(.est[.np])
  # delta method for the variance
  .seOut <- stats::setNames(rep(NA_real_, .np), .m$par)
  .seOut[.free] <- .se
  if (.m$par[.np] %in% .free) {
    .seOut[.np] <- .seOut[.np] * .est[.np]
  }
  list(est = .est, se = .seOut, m2ll = 2 * .o$value, convergence = .o$convergence)
}
