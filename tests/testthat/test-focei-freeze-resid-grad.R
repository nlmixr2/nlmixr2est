# freezeResidGrad: the outer FD gradient freezes the ODE (and the EBEs) when
# perturbing a residual/error theta.  These tests assert the option/accessor
# plumbing, that the freeze path actually fires (a silent fallback must fail
# loudly), and that a frozen fit matches the exact re-solve fit within tolerance.

nmTest({

  .freezeMod <- function() {
    ini({
      tka <- 0.45; tcl <- 1; tv <- 3.45
      eta.ka ~ 0.6; eta.cl ~ 0.3; eta.v ~ 0.1
      add.sd <- 0.7; prop.sd <- 0.1
    })
    model({
      ka <- exp(tka + eta.ka); cl <- exp(tcl + eta.cl); v <- exp(tv + eta.v)
      cp <- linCmt(); cp ~ add(add.sd) + prop(prop.sd)
    })
  }

  test_that("freezeResidGrad control default + residual-theta accessor", {
    # decided empirically (see NEWS): the freeze leaves the answer essentially
    # unchanged, so it ships on by default
    expect_true(foceiControl()$freezeResidGrad)
    expect_false(foceiControl(freezeResidGrad = FALSE)$freezeResidGrad)

    ui <- rxode2::rxode2(.freezeMod)
    # 0-based fullTheta indices of the non-fixed err thetas (add.sd=4, prop.sd=5)
    expect_equal(ui$foceiResidTheta, c(3L, 4L))
  })

  test_that("freeze fires for residual thetas and matches the exact gradient", {
    skip_on_cran()
    skip_if_not_installed("nlmixr2data")
    d <- nlmixr2data::theo_sd

    fit <- function(freeze) {
      f <- suppressMessages(suppressWarnings(
        nlmixr2(.freezeMod, d, "focei",
                foceiControl(freezeResidGrad = freeze, print = 0L,
                             covMethod = "", calcTables = FALSE))))
      list(ofv = as.numeric(f$objDf$OBJF), theta = fixef(f),
           nfreeze = nlmixr2est:::foceiNFreezeResidGrad())
    }

    off <- fit(FALSE)
    on  <- fit(TRUE)

    # mechanism actually ran (would be 0 on a silent fallback), and is genuinely
    # gated by the option
    expect_gt(on$nfreeze, 0)
    expect_equal(off$nfreeze, 0)

    # freezing is a small approximation, not bit-identical: objective and every
    # population estimate stay close to the exact re-solve
    expect_equal(on$ofv, off$ofv, tolerance = 1e-3)
    expect_equal(unname(on$theta), unname(off$theta), tolerance = 0.02)
  })

})
