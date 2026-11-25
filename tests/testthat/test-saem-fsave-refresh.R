nmTest({
  # The M-step resets phiM's phi0 columns to the new population value; the next
  # E-step scored the current state from predictions cached at the old value.
  .ctl <- saemControl(
    nBurn = 10,
    nEm = 5,
    seed = 42,
    print = 0,
    covMethod = "",
    calcTables = FALSE,
    # the non-centered steps move the chains themselves (and re-solve for it);
    # this file pins the phi0 broadcast alone
    muNonCentered = FALSE,
    omegaNonCentered = FALSE
  )

  test_that("SAEM refreshes cached predictions when a theta without an eta moves", {
    .mod <- function() {
      ini({
        tka <- 0.45
        tcl <- 1
        tv <- 3.45
        eta.cl ~ 0.3
        eta.v ~ 0.1
        add.sd <- 0.7
      })
      model({
        ka <- exp(tka)
        cl <- exp(tcl + eta.cl)
        v <- exp(tv + eta.v)
        linCmt() ~ add(add.sd)
      })
    }
    .n0 <- saemFsaveRefresh_()
    .fit <- .nlmixr(.mod, nlmixr2data::theo_sd, est = "saem", control = .ctl)
    .n1 <- saemFsaveRefresh_()
    expect_gt(.n1[["n"]], .n0[["n"]])
    # the cache it replaced was stale
    expect_gt(.n1[["drift"]], 0)
    expect_true(is.finite(.fit$objf))
  })

  test_that("SAEM does not re-solve when every theta has an eta", {
    .mod <- function() {
      ini({
        tka <- 0.45
        tcl <- 1
        tv <- 3.45
        eta.ka ~ 0.6
        eta.cl ~ 0.3
        eta.v ~ 0.1
        add.sd <- 0.7
      })
      model({
        ka <- exp(tka + eta.ka)
        cl <- exp(tcl + eta.cl)
        v <- exp(tv + eta.v)
        linCmt() ~ add(add.sd)
      })
    }
    .n0 <- saemFsaveRefresh_()
    .fit <- .nlmixr(.mod, nlmixr2data::theo_sd, est = "saem", control = .ctl)
    expect_equal(saemFsaveRefresh_()[["n"]], .n0[["n"]])
  })
})
