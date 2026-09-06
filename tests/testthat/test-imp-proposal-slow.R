# Fit-level coverage for impmapControl(proposal=).  The closed forms are pinned
# against quadrature in the essential test-imp-proposal.R; this file is the
# end-to-end half and lives in a weekly batch because it fits.
#
# The load-bearing test is the FIRST one: every family is an unbiased estimator
# of the SAME marginal likelihood, so if a family's draw and its density
# disagree the weights are wrong and the objective moves.  That is the only
# check that catches a draw/density mismatch -- a fit under a wrong density
# still converges to something plausible.
nmTest({
  .one <- function() {
    ini({
      tka <- 0.45; tcl <- 1; tv <- 3.45
      eta.ka ~ 0.6; eta.cl ~ 0.3
      add.sd <- 0.7
    })
    model({
      ka <- exp(tka + eta.ka)
      cl <- exp(tcl + eta.cl)
      v <- exp(tv)
      linCmt() ~ add(add.sd)
    })
  }

  test_that("every family estimates the same marginal likelihood", {
    .dat <- nlmixr2data::theo_sd
    # IS -2LL at the INITIAL parameters -- one iteration's first objective, so
    # every family is integrating the identical function.
    .obj1 <- function(seed, ...) {
      .f <- suppressWarnings(suppressMessages(
        nlmixr2(.one, .dat, "impmap",
                impmapControl(print = 0L, nIter = 1L, isample = 8000L,
                              impSeed = seed, auto = FALSE, covMethod = "",
                              calcTables = FALSE, ...))))
      .f$env$impObjTrace[1]
    }
    .seeds <- c(11L, 22L, 33L)
    .m <- c(
      normal  = mean(vapply(.seeds, function(s) .obj1(s), numeric(1))),
      t8      = mean(vapply(.seeds, function(s) .obj1(s, df = 8, proposal = "t"),
                            numeric(1))),
      laplace = mean(vapply(.seeds, function(s) .obj1(s, proposal = "laplace"),
                            numeric(1))),
      mixture = mean(vapply(.seeds, function(s) .obj1(s, proposal = "mixture"),
                            numeric(1)))
    )
    # A wrong density (e.g. a mixture weighted by the drawn component instead of
    # the mixture) biases this by order 0.2-4.6; Monte-Carlo noise at this
    # isample is order 0.01.  The bound is deliberately far below the former and
    # comfortably above the latter.
    expect_lt(diff(range(.m)), 0.15)
  })

  test_that("each family runs a full fit and reports itself", {
    .dat <- nlmixr2data::theo_sd
    .run <- function(...) {
      suppressWarnings(suppressMessages(
        nlmixr2(.one, .dat, "impmap",
                impmapControl(print = 0L, nIter = 6L, isample = 300L,
                              covMethod = "", calcTables = FALSE, ...))))
    }
    for (.p in c("normal", "laplace", "mixture")) {
      .f <- .run(proposal = .p)
      # the MECHANISM: the family the kernel actually used, per subject
      expect_identical(.f$env$impProposal, .p)
      if (identical(.p, "normal")) {
        # auto is on by default and the k-hat ladder lives on the normal/t
        # axis, so a subject may legitimately have been escalated -- that is
        # exactly what impPropInd exists to make visible.  Assert the only two
        # families reachable, and that any escalation agrees with impDfInd.
        expect_true(all(.f$env$impPropInd %in% c("normal", "t")))
        # impDfInd comes back as an nExp x 1 matrix (wrap of an arma::vec) and
        # impPropInd as a plain vector, so compare as vectors
        expect_identical(as.vector(.f$env$impPropInd == "t"),
                         as.vector(.f$env$impDfInd > 0))
      } else {
        # AUTO's df ladder is gated off for the non-df families, so no subject
        # may be converted out from under an explicit request
        expect_true(all(.f$env$impPropInd == .p))
        expect_true(all(.f$env$impDfInd == 0))
      }
      expect_true(is.finite(.f$objf))
      expect_true(all(is.finite(.f$env$impGammaInd)))
      # the sampler stayed usable
      expect_lt(max(.f$env$impPsisK), 0.7)
    }
    # the mixture reports its EFFECTIVE (covariance-matched) scales, which are
    # not the ones supplied
    .fm <- .run(proposal = "mixture", propMixScale = c(1, 9),
                propMixWeight = c(0.9, 0.1))
    expect_length(.fm$env$impPropMixScale, 2L)
    expect_equal(sum(.fm$env$impPropMixWeight * .fm$env$impPropMixScale), 1,
                 tolerance = 1e-10)
  })

  test_that("auto leaves a non-df family alone but still moves the budget", {
    .dat <- nlmixr2data::theo_sd
    .f <- suppressWarnings(suppressMessages(
      nlmixr2(.one, .dat, "impmap",
              impmapControl(print = 0L, nIter = 5L, isample = 300L,
                            proposal = "laplace", auto = TRUE,
                            covMethod = "", calcTables = FALSE))))
    # the df ladder is defined only on the normal/t axis, so no subject was
    # converted -- and auto=TRUE is NOT an error on a non-df family (it is the
    # default, so erroring would force every laplace user to pass auto=FALSE)
    expect_true(all(.f$env$impPropInd == "laplace"))
    expect_true(all(.f$env$impDfInd == 0))
    # the family-agnostic half of auto is still live
    expect_true(all(is.finite(.f$env$impNsampleInd)))
  })

  test_that("qrpem drives the new families too", {
    .dat <- nlmixr2data::theo_sd
    .f <- suppressWarnings(suppressMessages(
      nlmixr2(.one, .dat, "qrpem",
              qrpemControl(print = 0L, nIter = 4L, isample = 300L,
                           proposal = "laplace", qrScramble = "owen",
                           covMethod = "", calcTables = FALSE))))
    expect_identical(.f$env$impProposal, "laplace")
    expect_identical(.f$env$impQrScramble, "owen")
    expect_true(.f$env$impQr)
    expect_true(is.finite(.f$objf))
  })

  test_that("covMethod='imp' works under a new family", {
    .dat <- nlmixr2data::theo_sd
    .f <- suppressWarnings(suppressMessages(
      nlmixr2(.one, .dat, "impmap",
              impmapControl(print = 0L, nIter = 4L, isample = 300L,
                            proposal = "laplace", covMethod = "imp",
                            calcTables = FALSE))))
    # the covariance mirror must run under the same family, not fall back
    expect_true(is.finite(.f$objf))
    expect_true(all(is.finite(sqrt(diag(.f$cov)))))
  })
})
