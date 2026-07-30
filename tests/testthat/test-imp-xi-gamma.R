# Tests for the NONMEM-style per-subject xi statistic and the individual-gamma
# controller built on it.
#
# xi (NM7 Technical Guide eq. 1.73/1.75) is the mean importance weight normalized
# so the proposal kernel is 1 at its center:
#
#   xi_i = (1/r) sum_k [pi(eta_k) / pi(center)] / e(eta_k)
#
# It equals 1 when the proposal matches the posterior shape, is < 1 when the
# proposal is over-dispersed, and is > 1 when the posterior has heavier tails
# than the Gaussian proposal.  It is a DIFFERENT statistic from the Kish
# effective sample size (impNeffFrac) and the two are not comparable.
nmTest({

  .xiModel <- function() {
    ini({
      tka <- 0.45; tcl <- 1; tv <- 3.45
      eta.cl ~ 0.1
      add.sd <- 0.7
    })
    model({
      ka <- exp(tka)
      cl <- exp(tcl + eta.cl)
      v <- exp(tv)
      linCmt() ~ add(add.sd)
    })
  }

  test_that("xi: per-subject xi and its trace are exposed and well-formed", {
    .d <- nlmixr2data::theo_sd
    .f <- suppressWarnings(nlmixr2(.xiModel, .d, "impmap",
                                   impmapControl(print = 0L, nIter = 5L,
                                                 isample = 300L, covMethod = "")))
    .E <- .f$env
    # one xi per subject, all finite and strictly positive
    expect_false(is.null(.E$impXi))
    expect_length(.E$impXi, length(unique(.d$ID)))
    expect_true(all(is.finite(.E$impXi)))
    expect_true(all(.E$impXi > 0))
    # the trace has one entry per EM iteration actually run, alongside the
    # Kish-ESS trace it must not be confused with
    expect_length(.E$impXiTrace, .E$impIter)
    expect_length(.E$impNeffFrac, .E$impIter)
    expect_true(all(is.finite(.E$impXiTrace)))
    # xi and the Kish ESS fraction are different statistics: on a real fit they
    # must not coincide, which is exactly why $runInfo has to name which is which
    expect_false(isTRUE(all.equal(.E$impXiTrace, .E$impNeffFrac)))
  })

  test_that("xi decreases monotonically as the proposal is widened (gamma)", {
    # This is the property that makes xi targetable by the individual-gamma
    # controller: widening the proposal must lower xi, monotonically.  Compared
    # at iteration 1 across fits so every run sees identical starting parameters.
    .d <- nlmixr2data::theo_sd
    .xi1 <- vapply(c(1.0, 2.0, 4.0), function(g) {
      .f <- suppressWarnings(nlmixr2(.xiModel, .d, "impmap",
                                     impmapControl(print = 0L, nIter = 1L,
                                                   isample = 500L, gamma = g,
                                                   covMethod = "")))
      .f$env$impXiTrace[1]
    }, numeric(1))
    expect_true(all(is.finite(.xi1)))
    expect_true(all(diff(.xi1) < 0))
  })

  test_that("xi is near 1 for a well-matched proposal on a near-Gaussian model", {
    # gamma = 1 makes the proposal the Laplace approximation itself, so on a
    # model whose individual posterior is close to Gaussian xi should sit near 1.
    # Deliberately a loose band -- this pins the scale/orientation of the
    # statistic (a sign error or a missing normalizer would blow it out by
    # orders of magnitude), not its precise value.
    .d <- nlmixr2data::theo_sd
    .f <- suppressWarnings(nlmixr2(.xiModel, .d, "impmap",
                                   impmapControl(print = 0L, nIter = 1L,
                                                 isample = 500L, gamma = 1.0,
                                                 covMethod = "")))
    expect_true(.f$env$impXiTrace[1] > 0.5)
    expect_true(.f$env$impXiTrace[1] < 2.0)
  })

})
