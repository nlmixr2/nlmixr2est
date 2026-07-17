nmTest({

  test_that("issue #470: a theta-reset restart does not raise a factor-ID assertion", {

    # The fixed tv is deliberately far from the truth, so this fit never
    # converges (it drifts, triggers repeated "Theta reset (ETA drift)" and a
    # restart).  On the restart .foceiFitInternal() re-validates the previous
    # attempt's etaObf, whose ID column is a factor of the original subject IDs
    # -- which used to abort the whole fit with
    #   Assertion on 'fitEnv$etaObj$ID' failed: Must be of type 'integer', not 'factor'.
    # masking the real (non-convergence) reason.  The fit is still allowed to
    # fail, but never with that spurious assertion.
    mod <- function() {
      ini({
        tka <- 0.457920524576555
        tcl <- 1.01463922110082
        tv <- fix(7.42588307259747)
        add.sd <- c(0, 0.693895455003917)
        eta.ka ~ 0.406807327013617
        eta.cl ~ 0.0684420360035659
        eta.v ~ 0.019933306613011
      })
      model({
        ka <- exp(tka + eta.ka)
        cl <- exp(tcl + eta.cl)
        v <- exp(tv + eta.v)
        cp <- linCmt()
        cp ~ add(add.sd)
      })
    }

    .res <- tryCatch(
      .nlmixr(mod, theo_sd, "focei", foceiControl(print = 0L)),
      error = function(e) e)

    .msg <- if (inherits(.res, "condition")) conditionMessage(.res) else ""
    expect_false(grepl("not 'factor'", .msg, fixed = TRUE))
    expect_false(grepl("Must be of type 'integer'", .msg, fixed = TRUE))
  })

})
