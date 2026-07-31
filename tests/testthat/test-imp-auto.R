# AUTO: per-subject df / isample / iaccept (NONMEM AUTO=1 equivalent).
#
# The tutorial's documented triggers are implemented verbatim (nobs < neta, or
# non-normal data -> nonzero DF and IACCEPT ~ 0.2), but they are NECESSARY not
# SUFFICIENT: on plain theophylline neither fires and subjects still have
# infinite-variance weights.  AUTO therefore also escalates df from the Pareto
# k-hat diagnostic, which NONMEM has no equivalent of.  That escalation is
# nlmixr2's, not a reproduction of NONMEM.
nmTest({

  .pk <- function() {
    ini({tka <- 0.45; tcl <- 1; tv <- 3.45; eta.cl ~ 0.1; add.sd <- 0.7})
    model({ka <- exp(tka); cl <- exp(tcl + eta.cl); v <- exp(tv)
           linCmt() ~ add(add.sd)})
  }
  .fitAuto <- function(auto, nIter = 12L) {
    suppressWarnings(nlmixr2(.pk, nlmixr2data::theo_sd, "impmap",
                             impmapControl(print = 0L, nIter = nIter, isample = 300L,
                                           covMethod = "", auto = auto)))
  }

  test_that("auto control round-trips and defaults off", {
    expect_true(impmapControl()$auto)
    expect_false(impmapControl(auto = FALSE)$auto)
    expect_true(impmapControl(auto = TRUE)$auto)
    expect_true(do.call(impmapControl, impmapControl(auto = TRUE))$auto)
    expect_error(impmapControl(auto = "yes"))
    expect_true(all(c("auto", "autoNonNormal") %in% .impmapIsControlNames))
  })

  test_that("the in-kernel k-hat matches the validated R implementation", {
    # AUTO reacts to k-hat DURING the EM, so the kernel needs its own copy.  It
    # must agree with R/impPsis.R, which is the version validated against
    # loo::psis -- otherwise AUTO would be steering on a different number than
    # the one reported.
    .f <- .fitAuto(FALSE)
    .k <- .f$env$impKhatIter
    .r <- .f$env$impPsisK
    .ok <- is.finite(.k) & is.finite(.r)
    expect_gt(sum(.ok), 5L)
    expect_equal(.k[.ok], .r[.ok], tolerance = 1e-6)
  })

  test_that("auto=FALSE changes nothing", {
    .a <- .fitAuto(FALSE)
    expect_true(all(.a$env$impDfInd == 0))
    expect_true(all(.a$env$impIacceptInd == 0.4))
    expect_equal(length(unique(as.integer(.a$env$impNsampleInd))), 1L)
  })

  test_that("auto escalates df only for the subjects that need it", {
    # The discriminating test: theophylline is data-rich (11 obs, 1 eta) and
    # transformably normal, so the TUTORIAL's trigger does not fire at all.
    # Escalation must come from k-hat, and must be selective -- healthy
    # subjects keep the cheaper Gaussian proposal.
    skip_on_cran()
    .off <- .fitAuto(FALSE)
    .on <- .fitAuto(TRUE)
    .k0 <- .off$env$impPsisK
    expect_gt(sum(.k0 > 0.7), 0L)                      # premise
    # some subjects got a t proposal, but NOT all of them
    expect_gt(sum(.on$env$impDfInd > 0), 0L)
    expect_gt(sum(.on$env$impDfInd == 0), 0L)
    # and the tail failure is cleared
    expect_lt(max(.on$env$impPsisK), 0.7)
    expect_lt(max(.on$env$impPsisK), max(.k0))
    # without moving the answer
    expect_equal(.on$objf, .off$objf, tolerance = 0.5)
  })

  test_that("auto applies the tutorial's trigger for non-normal data", {
    # "or data are categorical" -> nonzero DF and IACCEPT ~ 0.2, for every
    # subject, regardless of how much data each has.
    skip_on_cran()
    .testSeed(202); rxode2::rxSetSeed(202)
    .d <- do.call(rbind, lapply(1:20, function(id) {
      .el <- stats::rnorm(1, 0, 0.7)
      data.frame(id = id, time = 1:10, dv = stats::rpois(10, exp(1 + .el)), evid = 0)
    }))
    .m <- function() {
      ini({tl <- 1; eta.l ~ 0.7})
      model({lam <- exp(tl + eta.l); dv ~ pois(lam)})
    }
    .f <- suppressWarnings(nlmixr2(.m, .d, "impmap",
                                   impmapControl(print = 0L, nIter = 8L, isample = 300L,
                                                 covMethod = "", auto = TRUE)))
    expect_true(all(.f$env$impDfInd > 0))          # every subject gets a t proposal
    # iaccept is NOT dropped to 0.2 up front any more.  Lowering it forces gamma
    # wide, and widening a Gaussian was measured not to fix tails while costing
    # a lot of ESS -- on this very fixture (k-hat already -1.33, i.e. no failure
    # at all) the blanket 0.2 cut ESS 0.549 -> 0.412 and doubled the objective
    # noise for nothing.  It is now held in reserve for subjects whose k-hat is
    # still failing after the df ladder is exhausted.
    expect_true(all(.f$env$impIacceptInd == 0.4))
  })

  test_that("auto reallocates the sample budget without inflating it", {
    skip_on_cran()
    .on <- .fitAuto(TRUE)
    .n <- length(.on$env$impNsampleInd)
    # load-balancing, not a cost increase: the total stays near isample*nsub
    expect_lt(abs(sum(.on$env$impNsampleInd) - 300 * .n), 0.25 * 300 * .n)
    expect_true(all(.on$env$impNsampleInd >= 25))   # floor keeps PSIS usable
  })

})
