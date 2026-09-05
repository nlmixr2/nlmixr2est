nmTest({
  # #1040: with mceta >= 1 the candidate set included the carried "last eta",
  # whose inner objective is essentially always the lowest, so mceta=n picked it
  # for every subject and returned the keep-last objective of mceta=-1/-2.  The
  # candidates are now eta=0 plus the (mceta-1) omega draws, and the converged
  # eta=0 solve is kept as a floor so mceta=n cannot end worse than mceta=0.
  #
  # The model needs curvature the inner problem can get lost in: two etas with a
  # FIXED unit omega entering through an inverse CDF, alongside an ordinary
  # estimated block.
  .mcetaMod <- function() {
    ini({
      tcl <- log(4)
      tv1 <- log(30)
      tq <- log(4)
      tv2 <- log(40)
      z.eta.cl ~ fix(1)
      z.eta.v1 ~ fix(1)
      eta.q + eta.v2 ~ c(0.0305,
                         0.0107, 0.0285)
      prop.sd <- 0.1
    })
    model({
      cl <- exp(tcl + 0.3 * logit(pnorm(z.eta.cl)))
      v1 <- exp(tv1 + 0.3 * logit(pnorm(z.eta.v1)))
      q <- exp(tq + eta.q)
      v2 <- exp(tv2 + eta.v2)
      linCmt() ~ prop(prop.sd)
    })
  }

  .mcetaData <- function() {
    withr::with_seed(42, {
      .ev <- rxode2::et(amt = 200, ii = 24, addl = 2)
      .ev <- rxode2::et(.ev, seq(0.5, 96, length.out = 12))
      .ev <- rxode2::et(.ev, id = 1:60)
      .d <- suppressWarnings(suppressMessages(
        as.data.frame(rxode2::rxSolve(.mcetaMod, .ev, addDosing = TRUE))))
    })
    .dat <- data.frame(ID = .d$id, TIME = .d$time, DV = .d$sim, EVID = .d$evid,
                       AMT = ifelse(is.na(.d$amt), 0, .d$amt))
    .dat$DV[.dat$EVID != 0] <- NA
    .dat
  }

  .mcetaFit <- function(dat, mceta, maxOuter = 0L, ...) {
    suppressWarnings(suppressMessages(
      nlmixr2(.mcetaMod, dat, "focei",
              foceiControl(print = 0L, covMethod = "", maxOuterIterations = maxOuter,
                           maxInnerIterations = 5000L, calcTables = FALSE,
                           mceta = mceta, ...))))
  }

  test_that("mceta > 0 explores its draws and cannot land worse than mceta=0", {
    skip_on_cran()
    .dat <- .mcetaData()
    .f0 <- .mcetaFit(.dat, 0L)
    .f1 <- .mcetaFit(.dat, 1L)
    .f5 <- .mcetaFit(.dat, 5L)

    # The counter is what shows the draws were USED.  Matching (or differing)
    # objectives cannot tell "explored and eta=0 won" from "never explored" --
    # the bug was exactly the second, and it looked like the first.
    expect_true(is.integer(.f5$env$nMcetaStart))
    expect_equal(sort(names(.f5$env$nMcetaStart)), c("sample", "zero"))
    expect_gt(.f5$env$nMcetaStart[["sample"]], 0L)

    # The floor comparison has to be made on the objective the fit REPORTS --
    # the marginal one, which carries the Laplace log|H| term -- not on the
    # inner joint density the optimizer minimizes.  The two order candidates
    # differently often enough that ranking on the inner objective alone hands
    # the fit the worse candidate: `flipped` counts the subjects where the
    # marginal disagreed, and it is a large fraction of `ranked` here.
    expect_gt(.f5$env$nInnerRerank[["ranked"]], 0L)
    expect_gt(.f5$env$nInnerRerank[["flipped"]], 0L)
    # A single candidate is not ranked at all, so mceta=0 pays nothing for this.
    expect_equal(.f0$env$nInnerRerank[["ranked"]], 0L)

    # Every inner problem has eta=0 among its candidates and keeps the converged
    # eta=0 solve as a floor, so at fixed parameters the objective cannot come
    # out above mceta=0.
    expect_lte(.f5$objf, .f0$objf + 1e-4)
    # ... and on this model it is strictly better, which is what "the extra
    # candidates are not being used" hid.
    expect_lt(.f5$objf, .f0$objf)

    # mceta=1 has no draws, so its candidate set is exactly {eta=0}.  It used to
    # be a silent no-op (empty sample cube -> the whole branch was skipped, so
    # the last eta was kept).
    expect_equal(.f1$objf, .f0$objf)
    expect_equal(.f1$env$nMcetaStart[["sample"]], 0L)
    expect_gt(.f1$env$nMcetaStart[["zero"]], 0L)
  })

  test_that("a fully mu-referenced model still honors mceta", {
    skip_on_cran()
    skip_if_not_installed("nlmixr2data")
    # A user-set mceta used to be silently reset to the default -2 for a model
    # whose etas are all mu-referenced, on the grounds that "the initial etas are
    # all exactly zero, so the search has nothing to explore".  That is only true
    # of the FIRST inner solve: every later one starts from the previous
    # iteration's mode, and the mceta>0 candidates are draws from omega, which
    # are not zero.  Resetting it is what made mceta=10 return an objective
    # bit-identical to mceta=-2 (#1040).
    .muMod <- function() {
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
    .warn <- character(0)
    .m5 <- withCallingHandlers(
      suppressMessages(
        nlmixr2(.muMod, nlmixr2data::theo_sd, "focei",
                foceiControl(print = 0L, covMethod = "", calcTables = FALSE,
                             mceta = 5L))),
      warning = function(w) {
        .warn <<- c(.warn, conditionMessage(w))
        invokeRestart("muffleWarning")
      })
    # The setting survived, and the draws were explored.
    expect_true(is.integer(.m5$env$nMcetaStart))
    expect_gt(.m5$env$nMcetaStart[["sample"]], 0L)
    # ... and nothing told the user the setting was being ignored.
    expect_false(any(grepl("has no effect", .warn, fixed = TRUE)))
  })

  test_that("the mceta search does not run on the outer gradient's difference legs", {
    skip_on_cran()
    # A finite-difference leg is pinned to the central evaluation's mode, so the
    # search has to be skipped there -- otherwise a tiny theta perturbation could
    # flip which candidate starts the inner problem and the two legs would sit in
    # different basins.  Count the inner solves that ran the search: one more
    # outer iteration costs a handful of objective evaluations (nsub starts
    # each), while an UNGATED search would add one full pass per parameter --
    # npars * nsub -- for that iteration's gradient alone.
    .dat <- .mcetaData()
    .nsub <- length(unique(.dat$ID))
    .g1 <- .mcetaFit(.dat, 5L, maxOuter = 1L, fast = FALSE)
    .g2 <- .mcetaFit(.dat, 5L, maxOuter = 2L, fast = FALSE)
    .npars <- length(fixef(.g1))
    expect_lt(sum(.g2$env$nMcetaStart) - sum(.g1$env$nMcetaStart), .nsub * .npars)
    # ... and the search did run somewhere, so the bound is not vacuous.
    expect_gt(.g1$env$nMcetaStart[["sample"]], 0L)
  })
})
