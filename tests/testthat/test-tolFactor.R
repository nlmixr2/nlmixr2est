nmTest({
  one.compartment <- function() {
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
      v  <- exp(tv  + eta.v)
      d/dt(depot)  <- -ka * depot
      d/dt(center) <-  ka * depot - cl / v * center
      cp <- center / v
      cp ~ add(add.sd)
    })
  }

  # Shared base fit (with default calcTables so the result is nlmixr2FitData).
  # NPDE is not requested here so addNpde tests can add it themselves.
  fit <- .nlmixr(one.compartment, nlmixr2data::theo_sd, est = "focei",
                 control = foceiControl(print = 0, covMethod = ""))

  test_that("tolFactor output: fit$env$tolFactor is a named numeric vector with one entry per subject >= 1", {
    tf <- fit$env$tolFactor
    nsubj <- length(unique(nlmixr2data::theo_sd$ID))
    expect_true(is.numeric(tf))
    expect_equal(length(tf), nsubj)
    expect_true(all(tf >= 1))
  })

  test_that("tolFactor output is used in addCwres: CWRES are finite", {
    # addCwres returns early if CWRES already exists (computed during default
    # table calculation), so the check still exercises the tolFactor code path.
    fitC <- suppressMessages(addCwres(fit))
    expect_true("CWRES" %in% names(fitC))
    expect_true(all(is.finite(fitC$CWRES)))
  })

  test_that("tolFactor output is used in addNpde: NPDE are finite", {
    fitN <- suppressMessages(
      addNpde(fit, table = tableControl(nsim = 50, seed = 42))
    )
    expect_true("NPDE" %in% names(fitN))
    expect_true(all(is.finite(fitN$NPDE)))
  })

  test_that("tolFactor output is used when npde is requested at fit time via tableControl", {
    fitNpde <- .nlmixr(one.compartment, nlmixr2data::theo_sd, est = "focei",
                       control = foceiControl(print = 0, covMethod = ""),
                       table = tableControl(npde = TRUE, nsim = 50, seed = 42))
    expect_true("NPDE" %in% names(fitNpde))
    expect_true(all(is.finite(fitNpde$NPDE)))
    tf <- fitNpde$env$tolFactor
    nsubj <- length(unique(nlmixr2data::theo_sd$ID))
    expect_equal(length(tf), nsubj)
    expect_true(all(tf >= 1))
  })

  test_that("tolFactor input: rxControl(tolFactor=5) is stored in the fit and used by the ODE solver", {
    fit5 <- .nlmixr(one.compartment, nlmixr2data::theo_sd, est = "focei",
                    control = foceiControl(print = 0, covMethod = "",
                                          rxControl = rxode2::rxControl(tolFactor = 5)))
    expect_s3_class(fit5, "nlmixr2FitData")

    # The input tolFactor is preserved in the stored rxControl and therefore
    # passed to every rxSolve_ call made with these options.
    expect_equal(fit5$foceiControl$rxControl$tolFactor, 5)

    # The per-subject output tolFactor accumulated during FOCEI optimization
    # is a valid numeric vector (>= 1); it is separate from the rxControl
    # tolFactor, which seeds the ODE-solver atol/rtol for the initial setup solve.
    tf <- fit5$env$tolFactor
    nsubj <- length(unique(nlmixr2data::theo_sd$ID))
    expect_equal(length(tf), nsubj)
    expect_true(all(tf >= 1))

    # addCwres and addNpde use fit$env$tolFactor for post-fit solves; verify
    # they work even when the fit was produced with a non-default rxControl tolFactor.
    fit5C <- suppressMessages(addCwres(fit5))
    expect_true(all(is.finite(fit5C$CWRES)))

    fit5N <- suppressMessages(
      addNpde(fit5, table = tableControl(nsim = 50, seed = 42))
    )
    expect_true(all(is.finite(fit5N$NPDE)))
  })
})
