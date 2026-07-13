# Mixture robustness guards (design/rpem/07, M2 gate): label-switching (components
# are exchangeable in the likelihood, so RPEM enforces an ascending-mu canonical
# order each iteration) and empty/collapsing components (a weight of exactly 0 makes
# logw = -Inf and permanently kills a component, so the M-step floors and renormalizes
# the weights).  These exercise K=2 mix() fits under conditions that stress both guards.

.genMix <- function(seed, tkaTrue, w1, nsub = 100L) {
  set.seed(seed); obsT <- seq(0.5, 24, by = 2)
  sim <- rxode2::rxode2({ ka <- exp(tka + eka); cl <- exp(tcl); v <- exp(tv); cp <- linCmt() })
  do.call(rbind, lapply(seq_len(nsub), function(i) {
    comp <- if (stats::runif(1) < w1) 1L else 2L
    ev <- rxode2::et(amt = 100, cmt = "depot"); ev <- rxode2::et(ev, obsT)
    s <- rxode2::rxSolve(sim, params = c(tka = tkaTrue[comp], tcl = 1, tv = 3.45,
                                         eka = stats::rnorm(1, 0, sqrt(0.3))),
                         events = ev, returnType = "data.frame", addDosing = FALSE)
    d <- as.data.frame(ev); d$id <- i; o <- d$evid == 0; d$DV <- 0
    d$DV[o] <- s$cp + stats::rnorm(nrow(s), 0, 0.1)
    d
  }))
}

.mixCtl <- function() rpemControl(nGauss = 500L, nMH = 80000L, mhBurn = 8000L,
                                  niter = 35L, collect = 12L, seed = 7L, cores = 4L)

test_that("RPEM mixture recovers unequal component weights (0.75/0.25)", {
  skip_on_cran()

  dat <- .genMix(42, c(0.0, 1.4), 0.75, nsub = 120L)
  mod <- function() {
    ini({ tka1 <- 0.2; tka2 <- 1.2; tcl <- fix(1.0); tv <- fix(3.45)
          p1 <- 0.5; add.sd <- 0.2; eta.ka ~ 0.3 })
    model({ ka <- mix(exp(tka1 + eta.ka), p1, exp(tka2 + eta.ka))
            cl <- exp(tcl); v <- exp(tv); cp <- linCmt(); cp ~ add(add.sd) })
  }
  rf <- .rpemFit(rxode2::rxode2(mod), dat, .mixCtl())
  expect_equal(unname(rf$mix$muK[1]), 0.0, tolerance = 0.3)
  expect_equal(unname(rf$mix$muK[2]), 1.4, tolerance = 0.3)
  # weights: wide bands (cores>1 MH is not bit-reproducible); w1+w2=1 so w[1] anchors it
  expect_true(rf$mix$w[1] > 0.6 && rf$mix$w[1] < 0.9)
  expect_equal(rf$addSd, 0.1, tolerance = 0.12)
})

test_that("RPEM mixture label-switching guard orders components (reversed init)", {
  skip_on_cran()

  dat <- .genMix(42, c(0.0, 1.4), 0.6, nsub = 100L)
  # start with the components in reversed order (tka1 > tka2); the ascending-mu guard
  # must still return components in canonical order matching the truth
  mod <- function() {
    ini({ tka1 <- 1.3; tka2 <- 0.1; tcl <- fix(1.0); tv <- fix(3.45)
          p1 <- 0.5; add.sd <- 0.2; eta.ka ~ 0.3 })
    model({ ka <- mix(exp(tka1 + eta.ka), p1, exp(tka2 + eta.ka))
            cl <- exp(tcl); v <- exp(tv); cp <- linCmt(); cp ~ add(add.sd) })
  }
  rf <- .rpemFit(rxode2::rxode2(mod), dat, .mixCtl())
  expect_lte(unname(rf$mix$muK[1]), unname(rf$mix$muK[2]))     # canonical (ascending) order
  expect_equal(unname(rf$mix$muK[1]), 0.0, tolerance = 0.25)
  expect_equal(unname(rf$mix$muK[2]), 1.4, tolerance = 0.25)
  expect_equal(unname(rf$mix$w[1]), 0.6, tolerance = 0.2)
})

test_that("RPEM mixture empty-component guard stays finite on a single population", {
  skip_on_cran()

  # data is a single population (both subjects drawn at tka = 0.7) but fit as K = 2;
  # the weight floor must keep every quantity finite (no logw = -Inf / NaN collapse)
  dat <- .genMix(11, c(0.7, 0.7), 1.0, nsub = 100L)
  mod <- function() {
    ini({ tka1 <- 0.3; tka2 <- 1.3; tcl <- fix(1.0); tv <- fix(3.45)
          p1 <- 0.5; add.sd <- 0.2; eta.ka ~ 0.3 })
    model({ ka <- mix(exp(tka1 + eta.ka), p1, exp(tka2 + eta.ka))
            cl <- exp(tcl); v <- exp(tv); cp <- linCmt(); cp ~ add(add.sd) })
  }
  rf <- .rpemFit(rxode2::rxode2(mod), dat, .mixCtl())
  expect_true(all(is.finite(c(rf$mix$muK, rf$mix$w, rf$omega, rf$addSd))))
  expect_true(all(rf$mix$w > 0))                              # no component died to zero
  expect_equal(sum(rf$mix$w), 1, tolerance = 1e-6)            # weights normalized
  expect_equal(rf$addSd, 0.1, tolerance = 0.12)               # residual still recovered
})
