# covMethod="analytic" with a record that is not an observation at the time of
# one: an EVID 2 record, or an EVID 0 record with a missing DV or MDV = 1 (both
# made EVID 2 by etTrans()).  The augmented sensitivity solves keep the solve's
# observation (EVID 0) records, so none of these may cost the analytic covariance.
#
# Two layers.  Full fits keep the analytic covariance and land near the fit without
# the record; the record adds a solver stop, which moves where the optimizer stops
# (largest measured: 2.4e-2 relative in a theta, 3.6e-2 in an SE, both for the
# record 0.05 h away).  The strict check evaluates the covariance at the base fit's
# estimates (maxOuterIterations = 0), where only the re-solved EBEs differ (largest
# measured SE difference: 3.8e-4 at the observation time, 9.8e-4 for the record
# 0.05 h away).  The subject solve itself, which compiles an augmented model, is the
# last test here.

nmTest({
  .orMod <- function() {
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
      d/dt(depot) <- -ka * depot
      d/dt(center) <- ka * depot - cl / v * center
      cp <- center / v
      cp ~ add(add.sd)
    })
  }

  # combined1 is not the plain add/prop case, so it takes the general (f,R) assembly
  .orModC1 <- function() {
    ini({
      tka <- 0.45
      tcl <- 1
      tv <- 3.45
      eta.ka ~ 0.6
      eta.cl ~ 0.3
      eta.v ~ 0.1
      add.sd <- 0.5
      prop.sd <- 0.1
    })
    model({
      ka <- exp(tka + eta.ka)
      cl <- exp(tcl + eta.cl)
      v <- exp(tv + eta.v)
      d/dt(depot) <- -ka * depot
      d/dt(center) <- ka * depot - cl / v * center
      cp <- center / v
      cp ~ add(add.sd) + prop(prop.sd) + combined1()
    })
  }

  .orBase <- nlmixr2data::theo_sd
  .orBase$MDV <- as.integer(.orBase$EVID != 0)
  # the third observation time of subject 1
  .orT <- .orBase$TIME[.orBase$ID == 1 & .orBase$EVID == 0][3]

  .orAdd <- function(time, evid, dv = NA_real_, mdv = 1L) {
    .r <- data.frame(
      ID = 1L,
      TIME = time,
      DV = dv,
      AMT = 0,
      EVID = evid,
      CMT = 2L,
      WT = .orBase$WT[.orBase$ID == 1][1],
      MDV = mdv
    )
    .d <- rbind(.orBase, .r)
    .d[order(.d$ID, .d$TIME, -.d$EVID), ]
  }

  .orCases <- list(
    evid2 = .orAdd(.orT, 2L),
    dvNa = .orAdd(.orT, 0L, NA_real_, 0L),
    mdv1 = .orAdd(.orT, 0L, 5, 1L),
    evid2Off = .orAdd(.orT + 0.05, 2L)
  )
  .orAtObs <- c("evid2", "dvNa", "mdv1")

  .orFit <- function(model, data, est = "focei", ...) {
    suppressMessages(suppressWarnings(nlmixr2(
      model,
      data,
      est,
      foceiControl(covMethod = "analytic", print = 0L, ...),
      table = tableControl(cwres = FALSE, npde = FALSE)
    )))
  }

  .orSe <- function(fit) sqrt(diag(fit$cov))

  .orExpectAnalytic <- function(fit) {
    expect_identical(fit$covMethod, "analytic (full)")
    expect_identical(fit$env$covTried, data.frame(method = "r", outcome = "used"))
  }

  # the covariance at `fit`'s estimates for each data set, against the base data's
  .orAtEstimates <- function(fit, est, tolAtObs, tolOff) {
    .b <- .orFit(fit$finalUi, .orBase, est, maxOuterIterations = 0L)
    .orExpectAnalytic(.b)
    .seB <- .orSe(.b)
    .ev <- lapply(.orCases, function(.d) .orFit(fit$finalUi, .d, est, maxOuterIterations = 0L))
    for (.n in names(.ev)) {
      .orExpectAnalytic(.ev[[.n]])
      expect_equal(.orSe(.ev[[.n]])[names(.seB)], .seB, tolerance = if (.n %in% .orAtObs) tolAtObs else tolOff)
    }
    # the three at-observation records are the same event once translated
    for (.n in .orAtObs[-1]) {
      expect_equal(.ev[[.n]]$objf, .ev$evid2$objf)
      expect_equal(.orSe(.ev[[.n]]), .orSe(.ev$evid2))
    }
    invisible(.ev)
  }

  test_that("a non-observation record at an observation time keeps the analytic covariance", {
    skip_on_cran()
    .base <- .orFit(.orMod, .orBase)
    .orExpectAnalytic(.base)
    # the first layer only says the refit stays in the same place: 0.05 in a theta
    # is twice the largest measured move, 0.1 in an SE nearly three times it
    for (.n in names(.orCases)) {
      .f <- .orFit(.orMod, .orCases[[.n]])
      .orExpectAnalytic(.f)
      expect_equal(.f$theta, .base$theta, tolerance = 0.05)
      expect_equal(.orSe(.f)[names(.orSe(.base))], .orSe(.base), tolerance = 0.1)
    }
    # at fixed estimates: 1e-3 and 2e-3 are 2.6x and 2x the largest measured differences
    .orAtEstimates(.base, "focei", tolAtObs = 1e-3, tolOff = 2e-3)
  })

  test_that("the general (f,R) assembly keeps it too, for FOCEI and FOCE", {
    skip_on_cran()
    for (.est in c("focei", "foce")) {
      .base <- .orFit(.orModC1, .orBase, .est)
      .orExpectAnalytic(.base)
      # largest measured: 3.8e-4 / 9.8e-4 (FOCEI) and 2.2e-4 / 5.1e-4 (FOCE)
      .orAtEstimates(.base, .est, tolAtObs = 1e-3, tolOff = 2e-3)
    }
  })

  test_that("fast = TRUE: the analytic gradient and the pooled covariance solve keep working", {
    skip_on_cran()
    .one <- function(data) {
      .n0 <- .odeSwapInfo()$pooledSolveN
      .b0 <- .foceiOuterFlagged$n
      .f <- .orFit(.orMod, data, fast = TRUE)
      list(fit = .f, pooled = .odeSwapInfo()$pooledSolveN - .n0, bailed = .foceiOuterFlagged$n - .b0)
    }
    .b <- .one(.orBase)
    .r <- .one(.orCases$evid2)
    for (.x in list(.b, .r)) {
      .orExpectAnalytic(.x$fit)
      expect_gt(.x$fit$env$nAnalyticGradDirect, 0)
      expect_gt(.x$pooled, 0L)
      expect_equal(.x$bailed, 0L)
    }
    # full refits, so the solver stop moves them (measured 8.9e-4 in a theta and
    # 9.3e-3 in an SE); 0.01 and 0.02 are about 11x and 2x that
    expect_equal(.r$fit$theta, .b$fit$theta, tolerance = 0.01)
    expect_equal(.orSe(.r$fit)[names(.orSe(.b$fit))], .orSe(.b$fit), tolerance = 0.02)
  })

  test_that("the subject solve keeps the observation records, whatever else shares their time", {
    skip_on_cran()
    .mod <- function() {
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
        d/dt(depot) <- -ka * depot
        d/dt(center) <- ka * depot - cl / v * center
        cp <- center / v
        cp ~ add(add.sd)
      })
    }
    .ui <- rxode2::rxUiDecompress(rxode2::assertRxUi(.mod))
    .am <- .foceiAnalyticAugModelDirs(.ui, paste0("ETA_", 1:3, "_"))
    expect_false(is.null(.am))
    .thr <- .ui$iniDf[!is.na(.ui$iniDf$ntheta), ]
    .thr <- .thr[order(.thr$ntheta), ]
    .p <- c(
      stats::setNames(.thr$est, paste0("THETA_", .thr$ntheta, "_")),
      stats::setNames(rep(0, 3), paste0("ETA_", 1:3, "_"))
    )
    .s1 <- nlmixr2data::theo_sd[nlmixr2data::theo_sd$ID == 1, ]
    .s1$MDV <- as.integer(.s1$EVID != 0)
    .t <- .s1$TIME[.s1$EVID == 0]
    .add <- function(evid, dv, mdv) {
      .r <- .s1[.s1$EVID == 0, ][1, ]
      .r$TIME <- .t[3]
      .r$EVID <- evid
      .r$DV <- dv
      .r$MDV <- mdv
      .d <- rbind(.s1, .r)
      .d[order(.d$TIME, -.d$EVID), ]
    }
    .e0 <- .foceiAnalyticSolveFA(.am, .p, .s1, .t)
    expect_length(.e0$f, length(.t))
    # the predictions agree to 1e-9: measured 0 for the missing-DV and MDV = 1
    # records and 1.1e-11 for the EVID 2 one, whose compartment event restarts the
    # solver at that time
    for (.d in list(.add(2L, NA_real_, 1L), .add(0L, NA_real_, 0L), .add(0L, 5, 1L))) {
      .e <- .foceiAnalyticSolveFA(.am, .p, .d, .t)
      expect_length(.e$f, length(.t))
      expect_equal(.e$f, .e0$f, tolerance = 1e-9)
    }
    # the times only confirm the rows; a mismatch is a NULL with the mismatch noted
    .foceiAnalyticWhyReset()
    expect_null(.foceiAnalyticSolveFA(.am, .p, .s1, c(.t, 30)))
    expect_identical(
      .foceiAnalyticWhy$detail,
      sprintf("%d observation rows solved for %d observations", length(.t), length(.t) + 1L)
    )
    # a failed rxSolve() is noted on one line
    .foceiAnalyticWhyReset()
    expect_null(.foceiAnalyticRxSolve(.am, .p[-1], .s1, 1e-8))
    expect_match(.foceiAnalyticWhy$detail, "^rxSolve: ")
    expect_false(grepl("\n", .foceiAnalyticWhy$detail, fixed = TRUE))
  })
})
