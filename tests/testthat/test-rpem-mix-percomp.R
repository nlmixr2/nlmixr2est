# Split-ETA per-component mixtures (design/rpem/07): each component carries its own
# random effect (eta.ka1 vs eta.ka2) and therefore its own Sigma^(k), declared as
# mix(exp(tka1 + eta.ka1), p1, exp(tka2 + eta.ka2)) with eta.ka1/eta.ka2 having their
# own omega.  RPEM draws each component's eta from its own omega, forces each component
# in the E-step, and estimates a per-component omega in the M-step.  On such data SAEM
# collapses the component means and reports only a single merged omega (its docs call
# split-component BSV "unreliable"); RPEM separates the components and reports one
# omega per component, so it meets/exceeds the SAEM mixture values here.

test_that("RPEM classifies a split-ETA per-component mixture", {
  mod <- function() {
    ini({ tka1 <- 0.2; tka2 <- 1.2; tcl <- fix(1.0); tv <- fix(3.45)
          p1 <- 0.5; add.sd <- 0.2; eta.ka1 ~ 0.2; eta.ka2 ~ 0.3 })
    model({ ka <- mix(exp(tka1 + eta.ka1), p1, exp(tka2 + eta.ka2))
            cl <- exp(tcl); v <- exp(tv); cp <- linCmt(); cp ~ add(add.sd) })
  }
  cl <- .rpemClassify(rxode2::rxode2(mod))
  expect_true(cl$mix$perComp)                              # each component has its own eta
  expect_equal(cl$mix$etaForComp, c(0L, 1L))              # component k -> eta k
  expect_equal(cl$nEta, 2L)
  expect_equal(cl$etaNames, c("eta.ka1", "eta.ka2"))
  expect_equal(unname(cl$muNames), c("tka1", "tka2"))
})

test_that("RPEM recovers a split-ETA per-component mixture with per-component Sigma", {
  skip_on_cran()
  skip_on_ci()  # heavy: multi-iteration mixture EM

  sim <- rxode2::rxode2({ ka <- exp(tka + eka); cl <- exp(tcl); v <- exp(tv); cp <- linCmt() })
  set.seed(42); nsub <- 150L; obsT <- seq(0.5, 24, by = 2)
  tkaTrue <- c(0.0, 1.4); omTrue <- c(0.1, 0.5); w1 <- 0.6
  dat <- do.call(rbind, lapply(seq_len(nsub), function(i) {
    k <- if (stats::runif(1) < w1) 1L else 2L
    ev <- rxode2::et(amt = 100, cmt = "depot"); ev <- rxode2::et(ev, obsT)
    s <- rxode2::rxSolve(sim, params = c(tka = tkaTrue[k], tcl = 1, tv = 3.45,
                                         eka = stats::rnorm(1, 0, sqrt(omTrue[k]))),
                         events = ev, returnType = "data.frame", addDosing = FALSE)
    d <- as.data.frame(ev); d$id <- i; o <- d$evid == 0; d$DV <- 0
    d$DV[o] <- s$cp + stats::rnorm(nrow(s), 0, 0.1)
    d
  }))
  mod <- function() {
    ini({ tka1 <- 0.2; tka2 <- 1.2; tcl <- fix(1.0); tv <- fix(3.45)
          p1 <- 0.5; add.sd <- 0.2; eta.ka1 ~ 0.2; eta.ka2 ~ 0.3 })
    model({ ka <- mix(exp(tka1 + eta.ka1), p1, exp(tka2 + eta.ka2))
            cl <- exp(tcl); v <- exp(tv); cp <- linCmt(); cp ~ add(add.sd) })
  }
  rf <- .rpemFit(rxode2::rxode2(mod), dat,
                 rpemControl(nGauss = 500L, nMH = 80000L, mhBurn = 8000L,
                             niter = 35L, collect = 12L, seed = 7L, cores = 4L))
  # components clearly separated (SAEM collapses this to ~0.45/0.67; the true gap is 1.4)
  expect_lt(rf$mix$muK[1], 0.5)
  expect_gt(rf$mix$muK[2], 1.0)
  expect_gt(rf$mix$muK[2] - rf$mix$muK[1], 1.0)
  # a distinct, finite, positive omega per component (SAEM reports only one merged omega)
  expect_equal(length(rf$omega), 2L)
  expect_true(all(is.finite(rf$omega) & rf$omega > 0))
  # weight and residual recovered
  expect_equal(unname(rf$mix$w[1]), 0.6, tolerance = 0.15)
  expect_equal(rf$addSd, 0.1, tolerance = 0.06)
})
