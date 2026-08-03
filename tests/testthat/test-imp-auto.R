# AUTO: per-subject df / isample / iaccept (NONMEM AUTO=1 equivalent).
#
# The tutorial's documented triggers are implemented verbatim (nobs < neta, or
# non-normal data -> nonzero DF and IACCEPT ~ 0.2), but they are NECESSARY not
# SUFFICIENT: on plain theophylline neither fires and subjects still have
# infinite-variance weights.  AUTO therefore also escalates df from the Pareto
# k-hat diagnostic, which NONMEM has no equivalent of.  That escalation is
# nlmixr2's, not a reproduction of NONMEM.
# NOTE ON gammaRule.  These tests exercise the TAIL machinery -- the t proposal
# (df), Pareto k-hat, and AUTO's df escalation -- all of which need a Gaussian
# proposal that actually FAILS in order to have anything to repair.  The default
# rule is now "target", which drives xi onto iaccept and in doing so repairs the
# tail itself: on the 3-ETA theophylline fixture it takes max k-hat from 0.836 to
# about -0.4, leaving these premises unsatisfiable.  So they pin
# gammaRule = "floor" deliberately.
#
# That overlap is a real consequence of the default change, not a test artifact:
# with "target" as the default the df/AUTO tail machinery is a secondary safety
# net rather than the primary remedy.
nmTest({

  # One eta on theophylline has NO tail failure (max k-hat about -1.8, nothing
  # above 0.7), so it cannot exercise AUTO's k-hat path at all.  What drives tail
  # failure is the number of ETAs, not the amount of data: identical data and
  # structural model with three etas reads max k-hat 1.13 with 3 of 12 subjects
  # failing.  .pk is kept for the "auto changes nothing" checks; .pk3 is the
  # fixture for anything asserting escalation.
  .pk <- function() {
    ini({tka <- 0.45; tcl <- 1; tv <- 3.45; eta.cl ~ 0.1; add.sd <- 0.7})
    model({ka <- exp(tka); cl <- exp(tcl + eta.cl); v <- exp(tv)
           linCmt() ~ add(add.sd)})
  }
  .pk3 <- function() {
    ini({tka <- 0.45; tcl <- 1; tv <- 3.45
         eta.ka ~ 0.3; eta.cl ~ 0.1; eta.v ~ 0.1; add.sd <- 0.7})
    model({ka <- exp(tka + eta.ka); cl <- exp(tcl + eta.cl); v <- exp(tv + eta.v)
           linCmt() ~ add(add.sd)})
  }
  # nobs < neta: the tutorial's documented sparse trigger.  Two observations per
  # subject against three etas, so the individual posterior is NOT identified and
  # the resulting heavy tail is structural -- no proposal shape repairs it.
  .sparseData <- local({
    set.seed(42)
    do.call(rbind, lapply(split(nlmixr2data::theo_sd, nlmixr2data::theo_sd$ID),
                          function(d) {
      .dose <- d[d$EVID != 0, , drop = FALSE]
      .obs <- d[d$EVID == 0, , drop = FALSE]
      rbind(.dose, .obs[sort(sample(seq_len(nrow(.obs)), 2L)), , drop = FALSE])
    }))
  })
  .fitAuto <- function(auto, nIter = 12L, model = .pk, data = nlmixr2data::theo_sd,
                       ...) {
    suppressWarnings(nlmixr2(model, data, "impmap",
                             impmapControl(print = 0L, nIter = nIter, isample = 300L,
                                           covMethod = "", auto = auto, ..., gammaRule = "floor")))
  }

  test_that("auto control round-trips and defaults off", {
    expect_true(impmapControl()$auto)
    expect_false(impmapControl(auto = FALSE, gammaRule = "floor")$auto)
    expect_true(impmapControl(auto = TRUE, gammaRule = "floor")$auto)
    expect_true(do.call(impmapControl, impmapControl(auto = TRUE, gammaRule = "floor"))$auto)
    expect_error(impmapControl(auto = "yes", gammaRule = "floor"))
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
    # The discriminating test.  Neither of the tutorial's triggers fires here
    # (11 observations against 3 etas, transformably normal), so escalation must
    # come from k-hat, and must be SELECTIVE -- healthy subjects keep the cheaper
    # Gaussian proposal.
    skip_on_cran()
    .off <- .fitAuto(FALSE, model = .pk3)
    .k0 <- .off$env$impPsisK
    # The premise, ASSERTED rather than skipped.  A skip would be quieter than
    # the thing it is guarding against: when the previous one-eta fixture stopped
    # producing tail failure, it was precisely this assertion going red that
    # surfaced it.  If this fails, the fixture no longer exercises AUTO and needs
    # re-selecting -- see plans/imp-auto-reinstrument.md -- rather than the
    # assertion being loosened.
    expect_gt(sum(.k0 > 0.7), 0L,
              label = "stressed-fixture subjects with k-hat > 0.7 (AUTO premise)")
    .on <- .fitAuto(TRUE, model = .pk3)
    # some subjects got a t proposal, but NOT all of them
    expect_gt(sum(.on$env$impDfInd > 0), 0L)
    expect_gt(sum(.on$env$impDfInd == 0), 0L)
    # and the tail failure is reduced
    expect_lt(max(.on$env$impPsisK), max(.k0))
    expect_lt(sum(.on$env$impPsisK > 0.7), sum(.k0 > 0.7))
    # without moving the answer
    expect_equal(.on$objf, .off$objf, tolerance = 0.5)
  })

  test_that("auto is close to free where there is no tail failure", {
    # The other half of "right defaults": on a model whose weights are already
    # well behaved, auto must cost essentially nothing.  This is the guard the
    # previous fixture's vacuous premise destroyed -- it went green while testing
    # nothing.
    #
    # Note what is NOT asserted: that nothing escalates.  impDfInd is the df of
    # the LAST E-step, but escalation is driven by the per-iteration k-hat, and
    # the early EM iterations run against a poorer proposal than the converged
    # one.  So a model that is healthy AT CONVERGENCE (max k-hat about -1.5,
    # nothing above 0.7) can still have escalated subjects on the transient.
    # Measured, that costs about 5% objective RMSE -- which is the claim worth
    # pinning, not a zero-escalation count that does not hold.
    skip_on_cran()
    .off <- .fitAuto(FALSE)
    .on <- .fitAuto(TRUE)
    expect_equal(sum(.off$env$impPsisK > 0.7), 0L)      # premise: converged fit healthy
    expect_gt(sum(.on$env$impDfInd == 0), 0L)           # escalation stays selective
    expect_equal(.on$objf, .off$objf, tolerance = 0.5)  # and does not move the answer
  })

  test_that("the nobs < neta trigger is gated, and autoNonmemSparse restores it", {
    # With fewer observations than random effects the individual posterior is not
    # identified, so the heavy tail is structural and no proposal shape repairs
    # it.  Applying the tutorial's rule there measurably hurts, so sparsity alone
    # must not assign a t proposal.
    skip_on_cran()
    .gated <- .fitAuto(TRUE, model = .pk3, data = .sparseData)
    .nonmem <- .fitAuto(TRUE, model = .pk3, data = .sparseData,
                        autoNonmemSparse = TRUE)
    # premise: every subject really is sparse in this fixture
    expect_true(all(.nonmem$env$impDfInd > 0))          # tutorial rule: everyone
    # gated: escalation is driven by k-hat, so it must not be universal-by-fiat
    expect_lt(sum(.gated$env$impDfInd > 0), sum(.nonmem$env$impDfInd > 0))
  })

  test_that("autoDfPatience controls withdrawal and round-trips", {
    expect_equal(impmapControl()$autoDfPatience, 2L)
    expect_equal(impmapControl(autoDfPatience = 0L, gammaRule = "floor")$autoDfPatience, 0L)
    expect_equal(do.call(impmapControl,
                         impmapControl(autoDfPatience = 3L, gammaRule = "floor"))$autoDfPatience, 3L)
    expect_error(impmapControl(autoDfPatience = -1L, gammaRule = "floor"))
    expect_error(impmapControl(autoDfPatience = "two", gammaRule = "floor"))
    expect_false(impmapControl()$autoNonmemSparse)
    expect_true(impmapControl(autoNonmemSparse = TRUE, gammaRule = "floor")$autoNonmemSparse)
    expect_error(impmapControl(autoNonmemSparse = "yes", gammaRule = "floor"))
    expect_true(all(c("autoNonmemSparse", "autoDfPatience") %in%
                      .impmapIsControlNames))
  })

  test_that("autoDfPatience = 0 keeps an escalation that patience would withdraw", {
    # patience 0 disables withdrawal.  On the structural fixture the shipped
    # default withdraws and 0 does not, so the two must differ -- otherwise the
    # control is inert and the withdrawal is not doing what it claims.
    skip_on_cran()
    .keep <- .fitAuto(TRUE, model = .pk3, data = .sparseData, autoDfPatience = 0L)
    .drop <- .fitAuto(TRUE, model = .pk3, data = .sparseData, autoDfPatience = 2L)
    expect_gte(sum(.keep$env$impDfInd > 0), sum(.drop$env$impDfInd > 0))
    expect_false(isTRUE(all.equal(.keep$objf, .drop$objf, tolerance = 1e-8)))
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
                                                 covMethod = "", auto = TRUE, gammaRule = "floor")))
    expect_true(all(.f$env$impDfInd > 0))          # every subject gets a t proposal
    # iaccept is NOT dropped to 0.2 up front any more.  Lowering it forces gamma
    # wide, and widening a Gaussian was measured not to fix tails while costing
    # a lot of ESS -- on this very fixture (k-hat already -1.33, i.e. no failure
    # at all) the blanket 0.2 cut ESS 0.549 -> 0.412 and doubled the objective
    # noise for nothing.  It is now held in reserve for subjects whose k-hat is
    # still failing after the df ladder is exhausted.
    expect_true(all(.f$env$impIacceptInd == 0.4))
    # Withdrawal must never take a non-normal subject BELOW the t proposal the
    # trigger assigned.  A non-normal endpoint needs heavy tails by construction,
    # and withdrawal is permanent, so a drop to Gaussian could not be undone --
    # the df floor is what prevents it.  Exercised at patience 1 to make
    # withdrawal as eager as it can be.
    .fp <- suppressWarnings(nlmixr2(.m, .d, "impmap",
                                    impmapControl(print = 0L, nIter = 8L, isample = 300L,
                                                  covMethod = "", auto = TRUE,
                                                  autoDfPatience = 1L, gammaRule = "floor")))
    expect_true(all(.fp$env$impDfInd > 0))
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
