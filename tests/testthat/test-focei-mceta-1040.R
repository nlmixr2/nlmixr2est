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
      rxz.eta.cl ~ fix(1)
      rxz.eta.v1 ~ fix(1)
      eta.q + eta.v2 ~ c(0.0305,
                         0.0107, 0.0285)
      prop.sd <- 0.1
    })
    model({
      cl <- exp(tcl + 0.3 * logit(pnorm(rxz.eta.cl)))
      v1 <- exp(tv1 + 0.3 * logit(pnorm(rxz.eta.v1)))
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

  .mcetaFit <- function(dat, mceta) {
    suppressWarnings(suppressMessages(
      nlmixr2(.mcetaMod, dat, "focei",
              foceiControl(print = 0L, covMethod = "", maxOuterIterations = 0L,
                           maxInnerIterations = 5000L, calcTables = FALSE,
                           mceta = mceta))))
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
})
