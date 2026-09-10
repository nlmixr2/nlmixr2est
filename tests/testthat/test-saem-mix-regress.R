test_that("saemControl(mixProbMethod='regress') fixes membership and separates the mixture", {
  skip_on_cran()

  # additive-error 2-component clearance mixture (5x separation), the same
  # well-behaved model used elsewhere in the mixture tests
  .testSeed(42)
  n_subj <- 30
  sub_pop <- rbinom(n_subj, 1, 0.6) + 1
  cl_sim <- ifelse(sub_pop == 1, 1.2, 6.0)
  sim_data <- do.call(rbind, lapply(seq_len(n_subj), function(i) {
    subj_cl <- cl_sim[i]; times <- c(0.5, 1, 2, 4, 8, 12, 24)
    ka_val <- 1.5; v_val <- 24.0; k_val <- subj_cl / v_val
    cp <- 100 * ka_val / (v_val * (ka_val - k_val)) *
      (exp(-k_val * times) - exp(-ka_val * times)) + rnorm(length(times), 0, 0.05)
    cp[cp < 0] <- 0
    data.frame(ID = i, TIME = c(0, times), AMT = c(100, rep(0, length(times))),
               EVID = c(1, rep(0, length(times))), DV = c(0, cp),
               CMT = c(1, rep(2, length(times))))
  }))

  mixmod <- function() {
    ini({
      tka <- log(1.5); tcl1 <- log(1.0); tcl2 <- log(5.0); tv <- log(20); p1 <- 0.5
      eta.cl ~ 0.01; eta.v ~ 0.01; eta.ka ~ 0.01; add.sd <- 0.05
    })
    model({
      ka <- exp(tka + eta.ka)
      cl <- mix(exp(tcl1 + eta.cl), p1, exp(tcl2 + eta.cl))
      v <- exp(tv + eta.v)
      d/dt(depot) <- -ka * depot
      d/dt(center) <- ka * depot - cl / v * center
      cp <- center / v
      cp ~ add(add.sd)
    })
  }

  fit <- suppressWarnings(nlmixr2(mixmod, sim_data, est = "saem",
    saemControl(print = 0, seed = 1234, nBurn = 250, nEm = 200,
                calcTables = FALSE, covMethod = 0L, mixProbMethod = "regress")))

  th <- fixef(fit)
  cls <- sort(exp(c(th[["tcl1"]], th[["tcl2"]])))
  # components must SEPARATE, not both collapse onto one value (the failure the
  # fixed-membership regressor prevents)
  expect_lt(cls[1], 3)
  expect_gt(cls[2], 3)
  # p1 is a finite proportion in (0,1)
  expect_true(is.finite(th[["p1"]]) && th[["p1"]] > 0 && th[["p1"]] < 1)
  # mechanism: fixed hard membership yields a per-subject classification
  expect_true(!is.null(fit$mixNum))
})

test_that("saem's mixture proportion equals its own responsibilities and the data's (#1058)", {
  skip_on_cran()

  # Two well-separated components (CL 1 vs 8, additive sd 0.05) so the realized
  # group fraction is the answer any method should find; everything except the
  # proportion and the BSV is fixed at truth.
  .testSeed(1001)
  nSub <- 60L
  clTrue <- c(1.0, 8.0)
  grp <- sample.int(2L, nSub, TRUE, prob = c(0.45, 0.55))
  sim <- rxode2::rxode2({
    ka <- 1.1
    cl <- CLI
    v <- 20
    d/dt(depot) <- -ka * depot
    d/dt(center) <- ka * depot - cl / v * center
    cp <- center / v
  })
  ev <- rxode2::et(rxode2::et(amt = 320, cmt = "depot"),
                   c(0.25, 0.5, 1, 2, 4, 8, 12, 24))
  obs <- do.call(rbind, lapply(seq_len(nSub), function(i) {
    s <- rxode2::rxSolve(sim, params = c(CLI = clTrue[grp[i]]), ev,
                         returnType = "data.frame")
    data.frame(ID = i, TIME = s$time, DV = s$cp + rnorm(nrow(s), 0, 0.05),
               AMT = 0, EVID = 0)
  }))
  dat <- rbind(data.frame(ID = seq_len(nSub), TIME = 0, DV = NA_real_,
                          AMT = 320, EVID = 1), obs)
  dat <- dat[order(dat$ID, dat$TIME, -dat$EVID), ]
  pTrue <- mean(grp == 1L)

  mixmod <- function() {
    ini({
      tka <- fix(log(1.1)); tcl1 <- fix(log(1.0)); tcl2 <- fix(log(8.0))
      tv <- fix(log(20)); p1 <- 0.45; eta.cl ~ 0.01; add.sd <- fix(0.05)
    })
    model({
      ka <- exp(tka)
      cl <- mix(exp(tcl1 + eta.cl), p1, exp(tcl2 + eta.cl))
      v <- exp(tv)
      linCmt() ~ add(add.sd)
    })
  }

  fit <- suppressWarnings(nlmixr2(mixmod, dat, est = "saem",
                                  saemControl(print = 0, nBurn = 100L, nEm = 100L,
                                              covMethod = "", calcTables = FALSE)))
  p1 <- fit$env$mixProbabilities[1]
  r1 <- fit$env$mixList[[1]]$prob

  # the proportion and the per-subject probabilities must describe ONE model:
  # the mixture score sum_i (r_i - p) is 0 at the EM/ML fixed point.  This is
  # data-free and is what caught the mlogit double back-transform.
  expect_equal(sum(r1 - p1), 0, tolerance = 1e-6)
  expect_equal(unname(p1), mean(r1), tolerance = 1e-8)
  expect_equal(sum(fit$env$mixProbabilities), 1, tolerance = 1e-8)
  # the reported theta is that same natural-scale proportion, not mexpit() of it
  expect_equal(unname(fit$theta[["p1"]]), unname(p1), tolerance = 1e-8)

  # ... and it is the proportion the data identifies (focei and imp both find
  # it on this dataset); saem reported expit(expit(p)) instead.
  expect_equal(unname(fit$env$mixProbabilities), c(pTrue, 1 - pTrue),
               tolerance = 0.05)
  expect_null(dim(fit$env$mixProbabilities))

  # mechanism: the proportion is the fraction of a per-subject classification,
  # so it is only right when that classification is.  Judging each component at
  # a subject's phiM draw -- whose fixed-effect-only columns carry a search
  # variance of 1, not a real BSV -- misclassified a fifth of these subjects
  # even though the components are 8-fold apart.
  expect_equal(as.integer(fit$env$mixNum$mixnum), as.integer(grp))

  # the identity has to survive the dedicated SA covariance phase too:
  # covMethod="sa" runs nSaCov extra iterations and then restores the converged
  # snapshot, and the reported proportion is taken after that restore
  fitSa <- suppressWarnings(nlmixr2(mixmod, dat, est = "saem",
                                    saemControl(print = 0, nBurn = 100L, nEm = 100L,
                                                covMethod = "sa", nSaCov = 50L,
                                                calcTables = FALSE)))
  .pSa <- fitSa$env$mixProbabilities[1]
  expect_equal(sum(fitSa$env$mixList[[1]]$prob - .pSa), 0, tolerance = 1e-6)
  expect_equal(unname(.pSa), pTrue, tolerance = 0.05)
})
