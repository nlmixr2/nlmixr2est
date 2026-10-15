# covMethod="analytic" with a record that is not an observation at the time of
# one: an EVID 2 record, or an EVID 0 record with a missing DV or MDV = 1 (both
# made EVID 2 by etTrans()).  The augmented sensitivity solves keep the solve's
# observation (EVID 0) records, so none of these may cost the analytic covariance.
#
# Two layers.  Full fits keep the analytic covariance and land near the fit without
# the record: an extra record adds a solver stop, which moves where the optimizer
# stops (measured up to 2.4e-2 relative in a theta and 3.6e-2 in an SE).  The strict
# check evaluates the covariance at the base fit's estimates (maxOuterIterations =
# 0), where only the re-solved EBEs differ (measured SEs within 4e-4 at the
# observation time, 1e-3 for the record 0.05 h away).

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
    for (.n in names(.orCases)) {
      .f <- .orFit(.orMod, .orCases[[.n]])
      .orExpectAnalytic(.f)
      expect_equal(.f$theta, .base$theta, tolerance = 0.05)
      expect_equal(.orSe(.f)[names(.orSe(.base))], .orSe(.base), tolerance = 0.05)
    }
    .orAtEstimates(.base, "focei", tolAtObs = 1e-3, tolOff = 2e-3)
  })

  test_that("the general (f,R) assembly keeps it too, for FOCEI and FOCE", {
    skip_on_cran()
    for (.est in c("focei", "foce")) {
      .base <- .orFit(.orModC1, .orBase, .est)
      .orExpectAnalytic(.base)
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
    expect_equal(.r$fit$theta, .b$fit$theta, tolerance = 0.01)
    expect_equal(.orSe(.r$fit)[names(.orSe(.b$fit))], .orSe(.b$fit), tolerance = 0.02)
  })

  test_that("the subject solve keeps the observation records, whatever else shares their time", {
    skip_on_cran()
    .ui <- rxode2::rxUiDecompress(rxode2::assertRxUi(.orMod))
    .am <- .foceiAnalyticAugModelDirs(.ui, paste0("ETA_", 1:3, "_"))
    expect_false(is.null(.am))
    .thr <- .ui$iniDf[!is.na(.ui$iniDf$ntheta), ]
    .thr <- .thr[order(.thr$ntheta), ]
    .p <- c(
      stats::setNames(.thr$est, paste0("THETA_", .thr$ntheta, "_")),
      stats::setNames(rep(0, 3), paste0("ETA_", 1:3, "_"))
    )
    .s1 <- function(d) d[d$ID == 1, , drop = FALSE]
    .t <- .orBase$TIME[.orBase$ID == 1 & .orBase$EVID == 0]
    .e0 <- .foceiAnalyticSolveFA(.am, .p, .s1(.orBase), .t)
    expect_length(.e0$f, length(.t))
    for (.n in .orAtObs) {
      .e <- .foceiAnalyticSolveFA(.am, .p, .s1(.orCases[[.n]]), .t)
      expect_length(.e$f, length(.t))
      expect_equal(.e$f, .e0$f, tolerance = 1e-6)
    }
    # the times only confirm the rows
    expect_null(.foceiAnalyticSolveFA(.am, .p, .s1(.orBase), c(.t, 30)))
  })
})
