## Psi conditional-likelihood builder (npBuildPsi / npEvalCondLik): validated
## against an INDEPENDENT rxode2 solve + Gaussian log-density.  The nonparametric
## weights are invariant to per-subject row scaling of Psi, so we compare log-Psi
## DIFFERENCES across support points per subject -- robust to the additive
## per-observation constant convention (e.g. the 2*pi normalizer).

nmTest({
  test_that("npBuildPsi matches an independent solve (log-Psi differences)", {
    npMod <- function() {
      ini({ tka <- log(1.5); tv <- log(32); tke <- log(0.08)
        eta.ka ~ 0.3; eta.v ~ 0.1; eta.ke ~ 0.1
        add.sd <- 0.7 })
      model({ ka <- exp(tka + eta.ka); v <- exp(tv + eta.v); ke <- exp(tke + eta.ke)
        d/dt(depot) <- -ka * depot
        d/dt(center) <- ka * depot - ke * center
        cp <- center / v
        cp ~ add(add.sd) })
    }
    ui <- rxode2::assertRxUi(npMod)
    ctl <- npagControl()
    dat <- nlmixr2data::theo_sd
    ids <- unique(dat$ID)
    N <- length(ids)

    # support points in eta space (eta.ka, eta.v, eta.ke)
    .testSeed(42)
    pts <- rbind(c(0, 0, 0),
                 c(0.4, -0.2, 0.1),
                 c(-0.5, 0.3, -0.15),
                 matrix(rnorm(9, 0, 0.3), 3, 3))
    colnames(pts) <- c("eta.ka", "eta.v", "eta.ke")

    .npInnerSetup(ui, dat, matrix(0, N, 3L), ctl)
    on.exit(.npInnerFree(), add = TRUE)
    psi <- .npInnerPsi(pts, ctl)
    expect_equal(dim(psi), c(N, nrow(pts)))
    expect_true(all(is.finite(psi)) && all(psi > 0))

    # Independent computation: solve the structural ODE at each support point and
    # sum the additive Gaussian log-densities (same -0.5 err^2/r - 0.5 log r form
    # that likInner0 accumulates; the 2*pi constant cancels in the differences).
    rxMod <- rxode2::rxode2({
      ka <- exp(tka + eta.ka); v <- exp(tv + eta.v); ke <- exp(tke + eta.ke)
      d/dt(depot) <- -ka * depot
      d/dt(center) <- ka * depot - ke * center
      cp <- center / v
    })
    th <- ui$theta
    addSd <- as.numeric(ui$iniDf$est[ui$iniDf$name == "add.sd"])
    rvar <- addSd * addSd

    indepLL <- matrix(NA_real_, N, nrow(pts))
    for (ii in seq_len(N)) {
      .d <- dat[dat$ID == ids[ii], , drop = FALSE]
      .obs <- .d[.d$EVID == 0, , drop = FALSE]
      for (kk in seq_len(nrow(pts))) {
        .p <- c(tka = unname(th[["tka"]]), tv = unname(th[["tv"]]), tke = unname(th[["tke"]]),
                eta.ka = unname(pts[kk, 1]), eta.v = unname(pts[kk, 2]), eta.ke = unname(pts[kk, 3]))
        .s <- rxode2::rxSolve(rxMod, params = .p, events = .d, returnType = "data.frame",
                              atol = 1e-10, rtol = 1e-10)
        # rxSolve returns one output row per observation record, in order
        .pred <- .s$cp
        expect_equal(length(.pred), nrow(.obs))
        .err <- .obs$DV - .pred
        indepLL[ii, kk] <- sum(-0.5 * .err * .err / rvar - 0.5 * log(rvar))
      }
    }

    # compare log-Psi differences (relative to support point 1) per subject.
    # log-Psi differences span O(10-100); the residual (~1e-2 absolute, ~1e-4
    # relative) is the inner solve's ODE tolerance, not a modeling difference.
    logPsi <- log(psi)
    dPsi <- logPsi - logPsi[, 1]
    dInd <- indepLL - indepLL[, 1]
    expect_lt(max(abs(dPsi - dInd)), 1e-2)
    # relative agreement across the dynamic range is tight
    expect_lt(max(abs(dPsi - dInd)) / max(abs(dInd)), 1e-3)
  })
})
