nmTest({
  test_that("mu-ref-covariate etas are protected from the eta-drift zero-reset", {
    # Phase 2 of the mu-referenced FOCEI family: mu-ref-covariate etas
    # (theta+eta+covariate, e.g. cl below) must never be zeroed by FOCEI's
    # internal eta-drift reset when muModel != "none", while plain etas
    # (theta+eta, no covariate, e.g. ka below) must still reset exactly as
    # today regardless of muModel.

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

    # Baseline (muModel="none"): the drift-reset fires normally for every
    # eta -- both columns collapse back toward zero after one inner step.
    expect_true(all(abs(etaNone$eta.cl) < 0.5))
    expect_true(all(abs(etaNone$eta.ka) < 0.5))

    # muModel="lin": the mu-ref-covariate eta (eta.cl) must be protected --
    # it stays much closer to its drifted starting point than the
    # unprotected baseline reset does. The plain eta (eta.ka) is NOT
    # mu-ref-covariate-eligible and must still reset like the baseline (a
    # single post-reset optimizer step can overshoot slightly past zero,
    # so compare against the baseline's own scale rather than a tight
    # absolute bound).
    expect_true(mean(abs(etaLin$eta.cl)) > 2 * mean(abs(etaNone$eta.cl)))
    expect_true(mean(abs(etaLin$eta.ka)) < 3 * mean(abs(etaNone$eta.ka)) + 0.2)
  })
})
