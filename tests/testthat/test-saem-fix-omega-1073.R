nmTest({
  # #1073: the reporting snapshot `Gamma2_phi1Report` used to be taken in the
  # middle of the omega M-step -- after the covstruct mask but BEFORE the
  # variance floors, the fix()ed-value restore and the diagonal enforcement.
  # The kernel sampled with the fixed value while `fit$omega` reported the
  # unconstrained M-step estimate.
  .mod <- function() {
    ini({
      tka <- 0.45
      tcl <- 1
      tv <- 3.45
      eta.ka ~ fix(0.3)
      eta.cl ~ 0.1
      add.sd <- 0.7
    })
    model({
      ka <- exp(tka + eta.ka)
      cl <- exp(tcl + eta.cl)
      v <- exp(tv)
      linCmt() ~ add(add.sd)
    })
  }

  .ctl <- saemControl(nBurn = 20, nEm = 20, print = 0, seed = 42,
                      calcTables = FALSE, covMethod = "")

  .fit <- suppressMessages(suppressWarnings(
    nlmixr2(.mod, nlmixr2data::theo_sd, "saem", .ctl)))

  test_that("saem reports a fix()ed eta variance as the fixed value (#1073)", {
    expect_equal(.fit$omega["eta.ka", "eta.ka"], 0.3, tolerance = 1e-12)
    # the unfixed variance is still estimated away from its ini()
    expect_true(.fit$omega["eta.cl", "eta.cl"] != 0.1)
  })

  test_that("the reporting snapshot is taken after the constraints (#1073)", {
    # mechanism: with a single component nothing is pooled, so the reported
    # matrix must be the live one the sampler used -- that equality is exactly
    # what the mis-ordered snapshot broke.
    expect_equal(.fit$saem$Gamma2_phi1Report, .fit$saem$Gamma2_phi1,
                 tolerance = 1e-12)
    expect_equal(.fit$saem$Gamma2_phi1Report[1, 1], 0.3, tolerance = 1e-12)
  })

  # A non-mu-referenced eta owns a phi column the same way, so it drifted the
  # same way.
  .modNonMu <- function() {
    ini({
      tka <- 0.45
      tcl <- 1
      tv <- 3.45
      eta.ka ~ fix(1)
      eta.cl ~ 0.1
      add.sd <- 0.7
    })
    model({
      ka <- exp(tka) * exp(eta.ka)
      cl <- exp(tcl + eta.cl)
      v <- exp(tv)
      linCmt() ~ add(add.sd)
    })
  }

  test_that("a non-mu-referenced fix()ed eta variance reports the fixed value (#1073)", {
    .fitNonMu <- suppressMessages(suppressWarnings(
      nlmixr2(.modNonMu, nlmixr2data::theo_sd, "saem", .ctl)))
    expect_equal(.fitNonMu$omega["eta.ka", "eta.ka"], 1, tolerance = 1e-12)
  })

  # A fix()ed CORRELATED block: the off-diagonal restore also happens after the
  # old snapshot point, so the reported covariance drifted too.
  .modBlock <- function() {
    ini({
      tka <- 0.45
      tcl <- 1
      tv <- 3.45
      add.sd <- 0.7
      eta.ka + eta.cl ~ fix(c(0.3, 0.02, 0.1))
    })
    model({
      ka <- exp(tka + eta.ka)
      cl <- exp(tcl + eta.cl)
      v <- exp(tv)
      linCmt() ~ add(add.sd)
    })
  }

  test_that("a fix()ed correlated omega block reports every fixed entry (#1073)", {
    .fitBlock <- suppressMessages(suppressWarnings(
      nlmixr2(.modBlock, nlmixr2data::theo_sd, "saem", .ctl)))
    expect_equal(unname(.fitBlock$omega),
                 matrix(c(0.3, 0.02, 0.02, 0.1), 2, 2),
                 tolerance = 1e-12)
  })

  # Reporting the constrained matrix means a floored variance now reaches
  # post-fit setup as a genuinely tiny (not negative) number, which can leave an
  # ESTIMATED off-diagonal at exactly 0.  rxSymInvCholCreate() reads that zero as
  # block structure and refuses the full cholesky parameter vector, so the omega
  # handed to it has to keep the declared pattern.
  .blockUi <- function() {
    ini({
      tka <- 0.45
      tcl <- 1
      tv <- 3.45
      add.sd <- 0.7
      eta.ka + eta.cl + eta.v ~ c(0.3, 0.02, 0.1, 0.01, 0.01, 0.1)
    })
    model({
      ka <- exp(tka + eta.ka)
      cl <- exp(tcl + eta.cl)
      v <- exp(tv + eta.v)
      linCmt() ~ add(add.sd)
    })
  }

  test_that("a declared omega off-diagonal that came back 0 stays representable (#1073)", {
    .ui <- rxode2::rxode2(.blockUi)
    .om <- .ui$omega
    .om[2, 3] <- .om[3, 2] <- 0
    expect_error(rxode2::rxSymInvCholCreate(mat = .om, diag.xform = "sqrt"))
    .fixed <- .foceiUnzeroDeclaredOffDiag(.om, .ui)
    expect_true(.fixed[2, 3] > 0)
    expect_equal(.fixed[2, 3], .fixed[3, 2])
    # every declared entry is back, so the cholesky vector matches again
    expect_equal(length(rxode2::rxSymInvCholCreate(mat = .fixed,
                                                   diag.xform = "sqrt")$theta), 6L)
    # still positive definite
    expect_false(inherits(try(chol(.fixed), silent = TRUE), "try-error"))
  })

  # A zero the model never declared is structure, not a degenerate estimate, and
  # must be left alone.
  .diagMod <- function() {
    ini({
      tka <- 0.45
      tcl <- 1
      tv <- 3.45
      add.sd <- 0.7
      eta.ka ~ 0.3
      eta.cl ~ 0.1
    })
    model({
      ka <- exp(tka + eta.ka)
      cl <- exp(tcl + eta.cl)
      v <- exp(tv)
      linCmt() ~ add(add.sd)
    })
  }

  test_that("a structural omega zero is not nudged (#1073)", {
    .diagUi <- rxode2::rxode2(.diagMod)
    .om <- .diagUi$omega
    expect_equal(.foceiUnzeroDeclaredOffDiag(.om, .diagUi), .om)
  })
})
