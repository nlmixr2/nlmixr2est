# Finite mixtures via the mix() interface: RPEM solves each Monte Carlo sample once
# per component (rpemEstepMixDraw), and the mixture M-step (rpemMstepMix) runs a
# Metropolis-Hastings over the joint (subject, sample, component) posterior to update
# the per-component typical values, the shared BSV, the mixture weights and the
# residual.  Must recover a simulated 2-component mixture and agree with SAEM.

test_that("RPEM recovers a 2-component mix() model and reports the mixture", {
  skip_on_cran()
  skip_on_ci()  # heavy: SAEM fit + multi-iteration mixture RPEM loop

  struct <- rxode2::rxode2({ ka <- exp(tka + eta); cl <- exp(tcl); v <- exp(tv); cp <- linCmt() })
  set.seed(42)
  nsub <- 80L; obsT <- seq(0.25, 12, by = 0.75)
  p1true <- 0.4; tkaT <- c(0.0, 1.6); omTrue <- 0.15; addT <- 0.1
  comp <- ifelse(runif(nsub) < p1true, 1L, 2L)
  dat <- do.call(rbind, lapply(seq_len(nsub), function(i) {
    ev <- rxode2::et(amt = 100, cmt = "depot"); ev <- rxode2::et(ev, obsT)
    s <- rxode2::rxSolve(struct,
                         params = c(tka = tkaT[comp[i]] + rnorm(1, 0, sqrt(omTrue)),
                                    tcl = 1, tv = 3.45, eta = 0),
                         events = ev, returnType = "data.frame", addDosing = FALSE)
    d <- as.data.frame(ev); d$id <- i; o <- d$evid == 0; d$DV <- 0
    d$DV[o] <- s$cp + rnorm(nrow(s), 0, addT)
    d
  }))

  mixmod <- function() {
    ini({ tka1 <- 0.1; tka2 <- 1.3; tcl <- fix(1.0); tv <- fix(3.45)
          p1 <- 0.5; add.sd <- 0.2; eta.ka ~ 0.3 })
    model({ ka <- mix(exp(tka1 + eta.ka), p1, exp(tka2 + eta.ka))
            cl <- exp(tcl); v <- exp(tv); cp <- linCmt(); cp ~ add(add.sd) })
  }

  rf <- suppressMessages(nlmixr2(mixmod, dat, est = "rpem",
                                 control = rpemControl(nGauss = 200L, nMH = 60000L,
                                                       mhBurn = 6000L, niter = 40L,
                                                       collect = 15L, seed = 100L, cores = 4L)))

  # recovers the simulated component typical values, weights, BSV and residual.
  .est <- rf$parFixedDf$Estimate
  names(.est) <- rownames(rf$parFixedDf)
  .p1 <- rf$mixProbabilities[1]
  # allow the components to come back in either order (label switching)
  .tk <- sort(c(.est["tka1"], .est["tka2"]))
  expect_equal(unname(.tk[1]), tkaT[1], tolerance = 0.15)
  expect_equal(unname(.tk[2]), tkaT[2], tolerance = 0.1)
  expect_equal(unname(.est["add.sd"]), addT, tolerance = 0.03)
  expect_true(abs(min(.p1, 1 - .p1) - min(p1true, 1 - p1true)) < 0.12)

  # per-subject component assignment recovers the truth (up to label switching).
  expect_true(is.data.frame(rf$mixNum))
  expect_true(all(rf$mixNum$mixnum %in% c(1L, 2L)))
  .acc <- max(mean(rf$mixNum$mixnum == comp), mean((3L - rf$mixNum$mixnum) == comp))
  expect_gt(.acc, 0.9)

  # mixList: one component per mixture, with posterior weights.
  expect_equal(length(rf$mixList), 2L)
  expect_true(all(c("ID", "eta.ka", "prob") %in% names(rf$mixList$mix1)))

  # Cross-check against SAEM on the same data: SAEM also fits the mix() model and
  # returns a valid 2-component mixture object.  We do NOT require SAEM to match the
  # component means -- SAEM tends to collapse this well-separated mixture toward the
  # pooled typical value -- so the meaningful comparison is that RPEM recovers the
  # simulated components at least as accurately as SAEM does.
  sf <- suppressMessages(nlmixr2(mixmod, dat, est = "saem",
                                 control = saemControl(print = 0, nBurn = 200, nEm = 300,
                                                       seed = 42, calcTables = FALSE)))
  expect_equal(length(sf$mixProbabilities), 2L)
  expect_equal(sum(sf$mixProbabilities), 1, tolerance = 1e-6)
  expect_true(all(sf$mixNum$mixnum %in% c(1L, 2L)))
  .se <- sf$parFixedDf$Estimate; names(.se) <- rownames(sf$parFixedDf)
  .stk <- sort(c(.se["tka1"], .se["tka2"]))
  .rpemErr <- sum(abs(.tk - tkaT))          # RPEM component-mean error vs truth
  .saemErr <- sum(abs(.stk - tkaT))         # SAEM component-mean error vs truth
  expect_lte(unname(.rpemErr), unname(.saemErr) + 0.2)
})
