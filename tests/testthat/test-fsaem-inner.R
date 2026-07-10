nmTest({
  # f-SAEM proposal core: reuse the FOCEi inner to get the per-subject conditional
  # MAP and the proposal precision H = Gamma_i^-1.  Validate that (a) the returned
  # eta is a stationary point of the individual objective and (b) H matches the
  # independent paper Eq-17 information J' J / sigma^2 + Omega^-1.
  test_that("fsaemInnerMap: MAP + Jacobian (Eq 17) proposal covariance", {
    skip_if_not_installed("numDeriv")
    skip_if_not_installed("nlmixr2data")

    one.cmt <- function() {
      ini({
        tka <- 0.45; tcl <- 1; tv <- 3.45
        eta.ka ~ 0.6; eta.cl ~ 0.3; eta.v ~ 0.1
        add.sd <- 0.7
      })
      model({
        ka <- exp(tka + eta.ka)
        cl <- exp(tcl + eta.cl)
        v  <- exp(tv + eta.v)
        linCmt() ~ add(add.sd)
      })
    }
    .data <- nlmixr2data::theo_sd
    .ui <- rxode2::rxUiDecompress(rxode2::rxode2(one.cmt))
    .neta <- 3L
    .N <- length(unique(.data$ID))

    .ctl <- list(rxControl = rxode2::rxControl(), fastCov = "jacobian", fastLik = "focei",
                 fastInnerIt = 100L, sumProd = FALSE, optExpression = TRUE, literalFix = FALSE,
                 addProp = "combined2", eventSens = "jump", indTolRelax = TRUE,
                 maxOdeRecalc = 5L, odeRecalcFactor = 10^0.5)

    .env <- .fsaemInnerSetup(.ui, .data, matrix(0, .N, .neta), .ctl)
    on.exit(.fsaemInnerFree(), add = TRUE)
    .map <- .fsaemInnerMap(.ctl, .neta)

    expect_equal(nrow(.map$eta), .N)
    expect_true(all(.map$ok == 1L))

    .sigma2 <- 0.7^2
    .OmegaInv <- solve(diag(c(0.6, 0.3, 0.1)))

    # prediction f_i(eta) via the inner driver (perturb only subject i's row)
    .predOf <- function(i, e, base) {
      .M <- base; .M[i, ] <- e
      as.numeric(vaeInnerLik(.M, 1L, FALSE, TRUE)$f[[i]])
    }
    # individual objective for subject i (for the MAP stationarity check)
    .objOf <- function(i, e) {
      .M <- .map$eta; .M[i, ] <- e
      vaeInnerLik(.M, 1L, FALSE, FALSE)$obj[i]
    }

    .maxRel <- 0
    .maxGrad <- 0
    for (i in seq_len(.N)) {
      .J <- numDeriv::jacobian(function(e) .predOf(i, e, .map$eta), .map$eta[i, ])
      .Hindep <- t(.J) %*% .J / .sigma2 + .OmegaInv
      .Hinner <- .map$hess[[i]]
      .maxRel <- max(.maxRel, max(abs(.Hinner - .Hindep) / pmax(abs(.Hindep), 1)))
      .maxGrad <- max(.maxGrad, max(abs(numDeriv::grad(function(e) .objOf(i, e), .map$eta[i, ]))))
    }
    # MAP is a stationary point of the individual objective
    expect_lt(.maxGrad, 1e-3)
    # proposal precision matches paper Eq 17 (relative; entries reach ~700)
    expect_lt(.maxRel, 1e-3)

    # proposal covariance is the inverse of the precision and is positive definite
    for (i in seq_len(.N)) {
      expect_true(all(is.finite(.map$gamma[[i]])))
      expect_true(all(eigen(.map$gamma[[i]], only.values = TRUE)$values > 0))
    }
  })
})
