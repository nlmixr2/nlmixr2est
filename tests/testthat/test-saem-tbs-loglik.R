nmTest({
  # saem's Gaussian-quadrature likelihood (calc.2LL) for a transform-both-sides
  # endpoint, cross-validated against a closed-form reference (#903).
  #
  # The reference is deliberately independent of saem: a 1-compartment IV bolus
  # has a closed-form prediction, and a single eta makes the marginal likelihood
  # a one-dimensional stats::integrate() written in plain R.  V is a model
  # constant rather than an estimated parameter so that every quantity the
  # quadrature uses comes straight from $theta/$omega, with no non-random phi
  # whose saem estimate can drift from the reported theta.

  .dose <- 320
  .v <- 70

  # log-normal residuals around dose/v*exp(-cl/v*t), one eta on CL
  .mkBolus <- function(seed, n, om, sd, times) {
    .testSeed(seed)
    .eta <- rnorm(n, 0, sqrt(om))
    do.call(rbind, lapply(seq_len(n), function(i) {
      .f <- .dose / .v * exp(-exp(log(4) + .eta[i]) / .v * times)
      rbind(data.frame(ID = i, TIME = 0, DV = NA_real_, AMT = .dose, EVID = 1),
            data.frame(ID = i, TIME = times,
                       DV = .f * exp(rnorm(length(times), 0, sd)),
                       AMT = 0, EVID = 0))
    }))
  }

  # additive residuals, for the untransformed control case
  .mkBolusAdd <- function(seed, n, om, sd, times) {
    .testSeed(seed)
    .eta <- rnorm(n, 0, sqrt(om))
    do.call(rbind, lapply(seq_len(n), function(i) {
      .f <- .dose / .v * exp(-exp(log(4) + .eta[i]) / .v * times)
      rbind(data.frame(ID = i, TIME = 0, DV = NA_real_, AMT = .dose, EVID = 1),
            data.frame(ID = i, TIME = times,
                       DV = .f + rnorm(length(times), 0, sd),
                       AMT = 0, EVID = 0))
    }))
  }

  # -2*log(marginal likelihood) of the ORIGINAL DV.  `tr` is the both-sides
  # transform and `ljac` its log-Jacobian log|dt/dy| -- what powerL supplies.
  # The defaults are lnorm; `ljac = .noJac` drops the Jacobian, giving the
  # likelihood on the transformed scale instead.
  #
  # The integrand is handled in the log domain, and the integration window is
  # placed around the mode using the local curvature: with many observations
  # the posterior is orders of magnitude narrower than the prior, and
  # integrating the raw density over the prior range makes stats::integrate()
  # miss the spike and return 0.  The per-subject modes come back in the
  # "logMax" attribute.
  .noJac <- function(y) rep(0, length(y))
  .refM2ll <- function(obs, tcl, om, sd, tr = log,
                       ljac = function(y) -log(y)) {
    .mx <- numeric(0)
    .ll <- vapply(unique(obs$ID), function(i) {
      .d <- obs[obs$ID == i, ]
      .lg <- function(e) {
        .f <- .dose / .v * exp(-exp(tcl + e) / .v * .d$TIME)
        sum(stats::dnorm(tr(.d$DV), tr(.f), sd, log = TRUE)) +
          sum(ljac(.d$DV)) + stats::dnorm(e, 0, sqrt(om), log = TRUE)
      }
      .lim <- 8 * sqrt(om)
      .grid <- seq(-.lim, .lim, length.out = 2001L)
      .k <- which.max(vapply(.grid, .lg, numeric(1)))
      .o <- stats::optimize(.lg,
                            c(.grid[max(1L, .k - 1L)],
                              .grid[min(length(.grid), .k + 1L)]),
                            maximum = TRUE, tol = 1e-12)
      .h <- 1e-5
      .s <- sqrt(-.h^2 / (.lg(.o$maximum + .h) - 2 * .o$objective +
                            .lg(.o$maximum - .h)))
      .mx <<- c(.mx, .o$objective)
      .o$objective + log(stats::integrate(function(eta) {
        exp(vapply(eta, .lg, numeric(1)) - .o$objective)
      }, max(-.lim, .o$maximum - 10 * .s), min(.lim, .o$maximum + 10 * .s),
      rel.tol = 1e-12, subdivisions = 2000L)$value)
    }, numeric(1))
    structure(-2 * sum(.ll), logMax = .mx)
  }

  .bolus <- function() {
    ini({
      tcl <- log(4)
      eta.cl ~ 0.09
      lnorm.sd <- 0.15
    })
    model({
      cl <- exp(tcl + eta.cl)
      v <- 70
      linCmt() ~ lnorm(lnorm.sd)
    })
  }

  .fitBolus <- function(d, nBurn = 40, nEm = 40) {
    suppressMessages(nlmixr2(.bolus, d, est = "saem",
      control = saemControl(nBurn = nBurn, nEm = nEm, seed = 1, print = 0L,
                            calcTables = FALSE)))
  }

  test_that(".logspaceAdd matches the naive sum and does not overflow", {
    .a <- c(-3, 0, 12.5, -Inf, -Inf, 800)
    .b <- c(2, 0, -4, 1.25, -Inf, 801)
    expect_equal(nlmixr2est:::.logspaceAdd(.a, .b),
                 c(log(exp(-3) + exp(2)), log(2), log(exp(12.5) + exp(-4)),
                   1.25, -Inf, log1p(exp(-1)) + 801))
    # the naive exp() form is Inf here; the log-domain form is finite
    expect_true(is.finite(nlmixr2est:::.logspaceAdd(800, 801)))
    expect_equal(nlmixr2est:::.logspaceAdd(800, 801),
                 nlmixr2est:::.logspaceAdd(801, 800))
  })

  test_that("saem calc.2LL for an lnorm() endpoint matches a closed-form reference (#903)", {
    .d <- .mkBolus(42L, 12L, 0.09, 0.15, c(0.5, 1, 2, 4, 7, 12, 24))
    .obs <- .d[.d$EVID == 0, ]
    .f <- .fitBolus(.d, nBurn = 60, nEm = 60)
    .ref <- as.numeric(.refM2ll(.obs, .f$theta[["tcl"]], .f$omega[1, 1],
                                .f$theta[["lnorm.sd"]]))
    .got <- suppressMessages(nlmixr2est:::calc.2LL(.f$saem, nnodes.gq = 25,
                                                  nsd.gq = 5, .f$phiM))
    # quadrature error is the only difference left (measured ~2e-3)
    expect_equal(.got, .ref, tolerance = 1e-4)

    # The defect was a SIGN: subtracting the transform's log-Jacobian instead of
    # adding it lands 4*sum(powerL) away, which for this data set is ~330 -- six
    # orders of magnitude outside the tolerance above.
    .wrong <- .ref - 4 * sum(log(.obs$DV))
    expect_gt(abs(.wrong - .ref), 100)

    # the transformed-scale likelihood is what the quadrature builds before the
    # Jacobian is applied; check the Jacobian moves it in the direction the
    # untransformed density requires
    .refNoJac <- as.numeric(.refM2ll(.obs, .f$theta[["tcl"]], .f$omega[1, 1],
                                     .f$theta[["lnorm.sd"]], ljac = .noJac))
    expect_equal(.ref, .refNoJac + 2 * sum(log(.obs$DV)), tolerance = 1e-4)

    # and the fit's own reported likelihood is the same corrected quantity (its
    # default quadrature is much coarser, hence the loose tolerance)
    expect_equal(-2 * as.numeric(logLik(.f)), .ref, tolerance = 2)
  })

  test_that("an untransformed add() endpoint is untouched by the Jacobian term (#903)", {
    # powerL is 0 for an untransformed endpoint, so the added term has to be
    # exactly zero and this likelihood has to be identical before and after the
    # sign fix.  Guards the other direction of #903: a Jacobian that leaks into
    # an add()/prop()/ll() fit.
    .addBolus <- function() {
      ini({
        tcl <- log(4)
        eta.cl ~ 0.09
        add.sd <- 0.2
      })
      model({
        cl <- exp(tcl + eta.cl)
        v <- 70
        linCmt() ~ add(add.sd)
      })
    }
    .d <- .mkBolusAdd(42L, 12L, 0.09, 0.2, c(0.5, 1, 2, 4, 7, 12, 24))
    .obs <- .d[.d$EVID == 0, ]
    .f <- suppressMessages(nlmixr2(.addBolus, .d, est = "saem",
      control = saemControl(nBurn = 60, nEm = 60, seed = 1, print = 0L,
                            calcTables = FALSE)))
    .ref <- as.numeric(.refM2ll(.obs, .f$theta[["tcl"]], .f$omega[1, 1],
                                .f$theta[["add.sd"]], tr = identity,
                                ljac = .noJac))
    .got <- suppressMessages(nlmixr2est:::calc.2LL(.f$saem, nnodes.gq = 25,
                                                  nsd.gq = 5, .f$phiM))
    expect_equal(.got, .ref, tolerance = 1e-5)
    # the transform is recorded as "no transform" (yj == 2), which is what makes
    # powerL return 0
    expect_equal(.f$saem$transMat[1, 2], 2)
  })

  test_that("the quadrature survives log-densities that overflow exp() (#903)", {
    # 500 observations per subject makes the per-subject log-density larger than
    # log(.Machine$double.xmax) = 709, so the old `Q <- Q + w*exp(ltot)`
    # accumulation overflowed and calc.2LL returned -Inf for the whole fit.
    .d <- .mkBolus(7L, 4L, 0.04, 0.08, seq(0.25, 40, length.out = 500L))
    .obs <- .d[.d$EVID == 0, ]
    .f <- .fitBolus(.d)
    .got <- suppressMessages(nlmixr2est:::calc.2LL(.f$saem, nnodes.gq = 9,
                                                  nsd.gq = 3, .f$phiM))
    .ref <- .refM2ll(.obs, .f$theta[["tcl"]], .f$omega[1, 1],
                     .f$theta[["lnorm.sd"]])
    # the exponent the accumulation has to carry (the mode of the log integrand,
    # plus the per-observation -0.5*log(2*pi) the quadrature leaves out) is well
    # past the overflow threshold
    expect_gt(max(attr(.ref, "logMax")) +
                0.5 * log(2 * pi) * max(table(.obs$ID)), 709)
    expect_true(is.finite(.got))
    expect_equal(.got, as.numeric(.ref), tolerance = 1e-4)
  })
})
