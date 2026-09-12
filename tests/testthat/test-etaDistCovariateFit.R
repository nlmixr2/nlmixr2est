# Phase 4.2 T5: the DEGENERATE covariate arm, end to end.
#
# Every other covariate test in test-etaDistCovariate.R is structural -- it
# asserts the parser sees the covariate, the decoder line carries it per record,
# the right thetas are held out.  None of them fits anything, so none of them
# would notice the estimator recovering a covariate effect that is not there.
#
# T5 simulates a covariate whose true coefficient is EXACTLY zero and requires
# the fit to say so.  That is the assertion a covariate search needs most: a
# method that reports an effect here would select spurious covariates on real
# data, and nothing structural can catch it.
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
    lclm <- 1.63; lv1m <- 1.55
    lclrv <- -2.4; lv1rv <- -2.4
    bWT <- 0
    eta.cl + eta.v1 ~ c(1, 0.5, 1)
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
  test_that("T5: a truly zero covariate effect is estimated as zero (focei)", {
    .d <- .edT5Data()
    .f <- suppressMessages(suppressWarnings(
      nlmixr2(.edT5Model(), .d, est = "focei",
              control = foceiControl(print = 0L, covMethod = ""))))
    .p <- setNames(.f$parFixedDf$Estimate, rownames(.f$parFixedDf))

    # the point of the arm: truth is exactly 0, and a real effect on this design
    # is 0.75 (T4), so 0.15 separates "found nothing" from "found something"
    expect_true(abs(.p[["bWT"]]) < 0.15)

    # prop.sd has lower=0 and a true value of 0.10.  It is asserted separately
    # because the recorded failure on an earlier version of this design was
    # prop.sd driven NEGATIVE (-0.889) and then a bounds violation -- a bounded
    # parameter leaving its bound is a defect whatever the covariate does.
    expect_true(.p[["prop.sd"]] > 0)
    expect_equal(.p[["prop.sd"]], 0.10, tolerance = 0.5)

    # and the rest of the model still lands on truth, so the zero above is a
    # real fit rather than an estimator that moved nothing
    expect_equal(.p[["lclm"]], 1.63, tolerance = 0.15)
    expect_equal(.p[["lv1m"]], 1.55, tolerance = 0.15)
    expect_true(.p[["lclrv"]] < -1.5)
  })
})
