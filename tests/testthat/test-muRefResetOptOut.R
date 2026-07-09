nmTest({
  test_that("mu-ref-covariate etas are protected from the eta-drift zero-reset", {
    # Phase 2 of the mu-referenced FOCEI family: mu-ref-covariate etas
    # (theta+eta+covariate, e.g. cl below) must be protected from FOCEI's internal
    # eta-drift reset when muModel != "none" (they are only ever updated by the
    # regression step), while plain etas (theta+eta, no covariate, e.g. ka below)
    # reset exactly as today regardless of muModel.

    theo_sd2 <- nlmixr2data::theo_sd
    theo_sd2$logWT <- log(theo_sd2$WT / 70)

    mod <- function() {
      ini({
        tka <- 0.45
        tcl <- 1
        tv <- 3.45
        allo.cl <- 0.75
        eta.ka ~ 0.6 # plain eta, always eligible for reset
        eta.cl ~ 0.3 # mu-ref-covariate eta, protected when muModel set
        eta.v ~ 0.1
        add.sd <- 0.7
      })
      model({
        ka <- exp(tka + eta.ka)
        cl <- exp(tcl + eta.cl + allo.cl * logWT)
        v <- exp(tv + eta.v)
        linCmt() ~ add(add.sd)
      })
    }

    # (1) WIRING: the reset-opt-out flag vector (foceiMuCovEtaVector, consumed by
    # src/inner.cpp isMuRefCovProtected) marks ONLY the mu-ref-covariate eta (eta.cl),
    # and only when muModel != "none".  etas are ordered by neta: ka, cl, v.
    .uiLin <- rxode2::rxUiDecompress(rxode2::rxode2(mod))
    rxode2::rxAssignControlValue(.uiLin, "muModel", "lin")
    expect_equal(unname(.uiLin$foceiMuCovEtaVector), c(0L, 1L, 0L))
    # muModel="none" (every non-mu method) sees the all-zero (empty) vector, unchanged
    .uiNone <- rxode2::rxUiDecompress(rxode2::rxode2(mod))
    expect_length(.uiNone$foceiMuCovEtaVector, 0L)

    # (2) RUNTIME: drift every eta far out, force the reset to fire (resetEtaP -> a tiny
    # standardized-eta bound), and take a single inner step.
    f0 <- .nlmixr(mod, theo_sd2, "focei", foceiControl(maxOuterIterations = 0L, print = 0))
    etaNames <- setdiff(names(f0$eta), "ID")
    etaMatDrift <- as.matrix(f0$eta[, etaNames, drop = FALSE])
    etaMatDrift[] <- 0
    etaMatDrift[, "eta.cl"] <- 2 # deliberately drifted (mu-ref-covariate)
    etaMatDrift[, "eta.ka"] <- 2 # deliberately drifted (plain)
    etaMatDrift[, "eta.v"] <- 0.01
    rownames(etaMatDrift) <- NULL

    runOne <- function(muModel) {
      fit <- .nlmixr(
        f0,
        est = "focei",
        control = foceiControl(
          maxOuterIterations = 0L, maxInnerIterations = 1L,
          etaMat = etaMatDrift, resetEtaP = 0.999,
          muModel = muModel, print = 0
        )
      )
      fit$eta
    }

    etaNone <- runOne("none")
    etaLin <- runOne("lin")

    # Baseline (muModel="none"): the drift-reset fires for every eta -- both columns
    # are pulled well back from their drifted start of 2.
    expect_true(mean(abs(etaNone$eta.cl)) < 1)
    expect_true(mean(abs(etaNone$eta.ka)) < 1)

    # muModel="lin": the mu-ref-covariate eta (eta.cl) is protected from the zero-reset,
    # so it retains proportionally MORE of its drift than the unprotected plain eta (ka)
    # does -- both measured against the muModel="none" baseline (the inner step converges
    # both toward the individual optimum, so the protection shows as a residual, not a
    # large absolute gap).
    .retainCl <- mean(abs(etaLin$eta.cl)) / mean(abs(etaNone$eta.cl))
    .retainKa <- mean(abs(etaLin$eta.ka)) / mean(abs(etaNone$eta.ka))
    expect_true(.retainCl > .retainKa + 0.1)
  })
})
