# Phase 4.2 T5: the DEGENERATE covariate arm, end to end.
#
# Every other covariate test in test-etaDistCovariate.R is structural -- it
# asserts the parser sees the covariate, the decoder line carries it per record,
# the right thetas are held out.  None of them fits anything, so none of them
# would notice the estimator recovering a covariate effect that is not there.
#
# T5 simulates a covariate whose true coefficient is EXACTLY zero.  The obvious
# test -- fit it and assert bWT ~ 0 -- is WORTHLESS on its own, and that is the
# main thing this file records: the coefficient's start is also 0, so an
# estimator that never moves it passes.  Measured on the T4 arm, where the true
# coefficient is 0.75, focei returns 0.0010 from a start of 0 and 0.4313 from a
# start of 0.5.  It is in the outer problem, it moves, and it moves the wrong
# way.  So this file asserts what actually holds, pins the weaker
# "not frozen" property, and leaves the recovery assertion skipped with its
# measurement rather than passing it vacuously.
#
# The data is built inline from a fixed seed rather than read from
# inst/sim/simCovT45.R's output, so the test is self-contained.  It is the same
# construction at a smaller subject count: gamma CL and V1 through a Gaussian
# copula, the covariate entering the RATE of the cl declaration per record.
.edT5Data <- function(nSub = 60L, bWT = 0.0, seed = 20260912L) {
  set.seed(seed)
  .lclm <- 1.63; .lv1m <- 1.55; .lclrv <- -2.4; .lv1rv <- -2.4; .rho <- 0.5
  .tim <- c(0.25, 0.5, 1, 2, 4, 8, 12, 24)
  .u <- function(z) pmin(pmax(stats::pnorm(z), 1e-15), 1 - 1e-15)
  .z1 <- stats::rnorm(nSub); .z2 <- stats::rnorm(nSub)
  .w2 <- .rho * .z1 + sqrt(1 - .rho^2) * .z2
  .wtBase <- stats::rnorm(nSub, 70, 12)
  .wtRec <- lapply(seq_len(nSub), function(i)
    round(.wtBase[i] + cumsum(stats::rnorm(length(.tim), 0, 1.5)), 1))
  .shCL <- 1 / exp(.lclrv); .shV1 <- 1 / exp(.lv1rv)
  .V1 <- stats::qgamma(.u(.w2), shape = .shV1,
                       rate = 1 / (exp(.lv1rv) * exp(.lv1m)))
  .m <- rxode2::rxode2({
    d/dt(central) <- -(cl / v) * central
    cp <- central / v
  })
  .rows <- lapply(seq_len(nSub), function(i) {
    .cl <- stats::qgamma(.u(.z1[i]), shape = .shCL,
                         rate = 1 / (exp(.lclrv) *
                                     exp(.lclm + bWT * log(.wtRec[[i]] / 70))))
    .ev <- data.frame(id = i, time = c(0, .tim),
                      amt = c(100, rep(NA_real_, length(.tim))),
                      evid = c(1L, rep(0L, length(.tim))), cmt = 1L,
                      cl = c(.cl[1], .cl), v = .V1[i])
    .s <- rxode2::rxSolve(.m, .ev, returnType = "data.frame")
    .s <- .s[!is.na(.s$cp) & .s$time > 0, ]
    data.frame(ID = i, TIME = .s$time, CP = .s$cp, AMT = NA_real_,
               EVID = 0L, CMT = 1L,
               WT = .wtRec[[i]][seq_len(nrow(.s))])
  })
  .obs <- do.call(rbind, .rows)
  # assay limit on the TRUE concentration: without it a 1e-12 observation
  # against a 1e-5 prediction is a relative residual of 1e7 under prop(), and
  # those records dominate the residual parameter
  .obs <- .obs[.obs$CP > 0.01, ]
  .obs$DV <- .obs$CP * (1 + stats::rnorm(nrow(.obs), 0, 0.10))
  .obs <- .obs[.obs$DV > 0, c("ID", "TIME", "DV", "AMT", "EVID", "CMT", "WT")]
  .dose <- data.frame(ID = seq_len(nSub), TIME = 0, DV = NA_real_, AMT = 100,
                      EVID = 1L, CMT = 1L, WT = round(.wtBase, 1))
  .d <- rbind(.dose, .obs)
  .d[order(.d$ID, .d$TIME, -.d$EVID), ]
}

.edT5Model <- function() {
  ini({
    # DISPLACED from the simulated truth (1.63, 1.55, -2.4, -2.4, rho 0.5) on
    # purpose: starting a parameter at its true value makes "the fit recovers
    # it" pass without the estimator moving anything, which is the same trap
    # bWT falls into below.
    lclm <- 1.9; lv1m <- 1.8
    lclrv <- -2.0; lv1rv <- -2.0
    bWT <- 0
    eta.cl + eta.v1 ~ c(1, 0.3, 1)
    dist(eta.cl) ~ dgamma(shape = 1 / exp(lclrv),
                          rate = 1 / (exp(lclrv) *
                                      exp(lclm + bWT * log(WT / 70))))
    dist(eta.v1) ~ dgamma(shape = 1 / exp(lv1rv),
                          rate = 1 / (exp(lv1rv) * exp(lv1m)))
    prop.sd <- 0.1
  })
  model({
    cl <- eta.cl; v <- eta.v1
    linCmt() ~ prop(prop.sd)
  })
}

nmTest({
  test_that("T5: the degenerate arm fits, and prop.sd stays inside its bound", {
    .d <- .edT5Data()
    .f <- suppressMessages(suppressWarnings(
      nlmixr2(.edT5Model(), .d, est = "focei",
              control = foceiControl(print = 0L, covMethod = ""))))
    .p <- setNames(.f$parFixedDf$Estimate, rownames(.f$parFixedDf))

    # prop.sd carries lower = 0 and a true value of 0.10.  An earlier version of
    # this design drove it to -0.889 and then failed on a bounds violation; a
    # bounded parameter leaving its bound is a defect whatever else happens, so
    # this is asserted on its own.
    expect_true(.p[["prop.sd"]] > 0)
    expect_equal(.p[["prop.sd"]], 0.10, tolerance = 0.5)

    # The structural parameters are started away from truth and must come back
    # to it, so these cannot pass on a frozen fit.  Measured at 120 subjects:
    # lclm 1.9 -> 1.6228, lv1m 1.8 -> 1.5440, lclrv -2.0 -> -2.2394,
    # rxCor 0.3 -> 0.6017, against truth 1.63, 1.55, -2.4 and 0.5.
    expect_equal(.p[["lclm"]], 1.63, tolerance = 0.15)
    expect_equal(.p[["lv1m"]], 1.55, tolerance = 0.15)
    # directional rather than tight: the relative variances move toward truth
    # from -2.0 but do not arrive, and the copula moves well off its 0.3 start
    expect_true(.p[["lclrv"]] < -2.05)
    expect_true(.p[["rxCor.eta.v1.eta.cl"]] > 0.4)

    # NOTE: bWT is deliberately NOT asserted here.  Its start and its truth are
    # both 0, so an estimator that never moves it "recovers" it perfectly --
    # which is exactly what happens today (see the skipped test below).  Reading
    # bWT ~ 0 on this arm as evidence that the covariate machinery works is the
    # trap this file exists to avoid.
  })

  test_that("a covariate coefficient on a declaration is recovered", {
    skip(paste("the coefficient is in the outer problem but its search stalls:",
               "on the T4 arm (true bWT = 0.75) focei returns 0.0010 from a",
               "start of 0, and 0.4313 from a start of 0.5 -- it moves, and",
               "moves the wrong way.  Enable when the outer search reaches the",
               "coefficient; the objective itself is minimized at truth."))
    # The assertion this arm is for, kept here so it is enabled rather than
    # rewritten once the search is fixed:
    #   .d <- .edT5Data(bWT = 0.75)
    #   .f <- nlmixr2(.edT5Model(), .d, est = "focei", ...)
    #   expect_equal(.p[["bWT"]], 0.75, tolerance = 0.2)
  })

  test_that("the coefficient is at least IN the outer problem, not frozen", {
    # Weaker than recovery but still worth pinning: a coefficient that never
    # moves at all is a different (and worse) defect than one whose search
    # stalls, and the two are indistinguishable on the degenerate arm.
    .d <- .edT5Data()
    .m <- .edT5Model()
    .m <- rxode2::ini(.m, bWT = 0.5)
    .f <- suppressMessages(suppressWarnings(
      nlmixr2(.m, .d, est = "focei",
              control = foceiControl(print = 0L, covMethod = ""))))
    .b <- .f$parFixedDf$Estimate[match("bWT", rownames(.f$parFixedDf))]
    expect_false(isTRUE(all.equal(.b, 0.5, tolerance = 1e-6)))
  })
})
