nmTest({
  test_that("output/posthoc declare full prior support (#938)", {
    # "output" and "posthoc" do not estimate anything; refusing a prior there
    # would only break assembling a finished fit from a prior-carrying model.
    for (.cls in c("output", "posthoc")) {
      .env <- new.env(parent = emptyenv())
      class(.env) <- c(.cls, "nlmixr2Est")
      expect_identical(.nlmixr2PriorSupport(.env), "all")
    }
    # a method without the attribute still defaults to "none" (focei's
    # family declares "theta" as of #931, so it is no longer a
    # representative "undeclared" example -- saem is)
    .env <- new.env(parent = emptyenv())
    class(.env) <- c("saem", "nlmixr2Est")
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
        d/dt(depot) <- -ka * depot
        d/dt(center) <- ka * depot - cl / v * center
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

    # est="focei" now declares nlmixr2Priors = "general" (#931): it honours
    # a prior on a population parameter AND on an omega element, so this is
    # no longer a case that demonstrates "the bypass is scoped" -- saem
    # (which declares no prior support at all) is.
    saem.prior <- function() {
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
        d/dt(depot) <- -ka * depot
        d/dt(center) <- ka * depot - cl / v * center
        cp <- center / v
        cp ~ add(add.sd)
      })
    }
    # the bypass is scoped: the gate still refuses a user-initiated saem
    expect_error(
      suppressWarnings(suppressMessages(
        nlmixr2(saem.prior, nlmixr2data::theo_sd, est = "saem"))),
      "prior")
    expect_false(isTRUE(nlmixr2global$nlmixr2PriorGateBypass))

    # posthoc/output declare nlmixr2Priors = "all" (#938); est="focei" now
    # also fully supports an omega-touching prior (#931) -- both correctly
    # evaluate one (a previous version of this fix left a segfault
    # reachable here: op_focei.omega is populated only for est="fo", so
    # reading it unconditionally read a default-constructed (0x0, NULL
    # memptr()) matrix for every other method; foceiCurrentOmega()
    # recovers a live Omega from op_focei.omegaInv instead).
    one.compartment.omega.prior <- function() {
      ini({
        tka <- 0.45
        tcl <- 1
        tv <- 3.45
        add.sd <- c(0, 0.7)
        eta.ka ~ 0.6
        prior(eta.ka) ~ dnorm(0, 0.3)
      })
      model({
        ka <- exp(tka + eta.ka)
        cl <- exp(tcl)
        v <- exp(tv)
        d/dt(depot) <- -ka * depot
        d/dt(center) <- ka * depot - cl / v * center
        cp <- center / v
        cp ~ add(add.sd)
      })
    }
    .fitOmega <- suppressWarnings(suppressMessages(
      nlmixr2(one.compartment.omega.prior, nlmixr2data::theo_sd, est = "posthoc")))
    expect_true(inherits(.fitOmega, "nlmixr2FitData"))
    expect_true("eta.ka" %in% rxode2::rxUiPriors(.fitOmega$ui)$name)
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
        d/dt(depot) <- -ka * depot
        d/dt(center) <- ka * depot - cl / v * center
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
