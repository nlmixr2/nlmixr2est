nmTest({
  test_that("output/posthoc declare full prior support (#938)", {
    # "output" and "posthoc" do not estimate anything; refusing a prior there
    # would only break assembling a finished fit from a prior-carrying model.
    for (.cls in c("output", "posthoc")) {
      .env <- new.env(parent = emptyenv())
      class(.env) <- c(.cls, "nlmixr2Est")
      expect_identical(.nlmixr2PriorSupport(.env), "all")
    }
    # a method without the attribute still defaults to "none"
    .env <- new.env(parent = emptyenv())
    class(.env) <- c("focei", "nlmixr2Est")
    expect_identical(.nlmixr2PriorSupport(.env), "none")
  })

  test_that("posthoc runs on a model whose ini({}) declares priors (#938)", {
    skip_on_cran()
    skip_if_not(exists("rxUiPriors", envir = asNamespace("rxode2")),
                "rxode2 without prior support")

    one.compartment <- function() {
      ini({
        tka <- 0.45
        tcl <- 1
        tv <- 3.45
        add.sd <- c(0, 0.7)
        eta.ka ~ 0.6
        prior(tka) ~ dnorm(0, 10)
        prior(add.sd) ~ dcauchy(0, 5)
      })
      model({
        ka <- exp(tka + eta.ka)
        cl <- exp(tcl)
        v <- exp(tv)
        d / dt(depot) <- -ka * depot
        d / dt(center) <- ka * depot - cl / v * center
        cp <- center / v
        cp ~ add(add.sd)
      })
    }

    # Before #938 this errored at the prior gate ("...which this estimation
    # method cannot use"): first at the posthoc dispatch itself, then again in
    # the automatic FOCEi-objective evaluation (.setOfvFo), which re-enters
    # with est="focei".
    .fit <- suppressWarnings(suppressMessages(
      nlmixr2(one.compartment, nlmixr2data::theo_sd, est = "posthoc")))
    expect_true(inherits(.fit, "nlmixr2FitData"))
    # the priors survive onto the finished fit
    .pri <- rxode2::rxUiPriors(.fit$ui)
    expect_true(all(c("tka", "add.sd") %in% .pri$name))
    # the bypass is scoped: the gate still refuses a user-initiated focei
    expect_error(
      suppressWarnings(suppressMessages(
        nlmixr2(one.compartment, nlmixr2data::theo_sd, est = "focei",
                control = foceiControl(maxOuterIterations = 0L, print = 0L)))),
      "prior")
    expect_false(isTRUE(nlmixr2global$nlmixr2PriorGateBypass))
  })

  test_that("addCwres() works on a prior-carrying fit (#938)", {
    skip_on_cran()
    skip_if_not(exists("rxUiPriors", envir = asNamespace("rxode2")),
                "rxode2 without prior support")

    one.compartment <- function() {
      ini({
        tka <- 0.45
        tcl <- 1
        tv <- 3.45
        add.sd <- c(0, 0.7)
        eta.ka ~ 0.6
        prior(tka) ~ dnorm(0, 10)
      })
      model({
        ka <- exp(tka + eta.ka)
        cl <- exp(tcl)
        v <- exp(tv)
        d / dt(depot) <- -ka * depot
        d / dt(center) <- ka * depot - cl / v * center
        cp <- center / v
        cp ~ add(add.sd)
      })
    }
    # cwres=FALSE keeps CWRES out of the fit so addCwres() genuinely runs its
    # nested est="focei" evaluation (with CWRES already present it returns
    # early and the bypass is never exercised)
    .fit <- suppressWarnings(suppressMessages(
      nlmixr2(one.compartment, nlmixr2data::theo_sd, est = "posthoc",
              table = tableControl(cwres = FALSE))))
    expect_false("CWRES" %in% names(.fit))
    .fit2 <- suppressWarnings(suppressMessages(addCwres(.fit)))
    expect_true("CWRES" %in% names(.fit2))
    expect_false(isTRUE(nlmixr2global$nlmixr2PriorGateBypass))
  })
})
