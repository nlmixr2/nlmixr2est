# Multiple mix() calls (design/rpem/07): more than one mixed parameter sharing ONE
# latent class -- the same number of components and the same mixture probabilities
# (e.g. ka and cl both mixed, one subpopulation label governing both).  RPEM parses
# each mix() call as a mixed parameter, forces each component in the E-step, and the
# M-step accumulates every parameter's per-component typical value under the shared
# class label.

test_that("RPEM classifies multiple mix() calls (shared latent class)", {
  mod <- function() {
    ini({ tka1 <- 0.2; tka2 <- 1.2; tcl1 <- 0.8; tcl2 <- 1.5; tv <- fix(3.45)
          p1 <- 0.5; add.sd <- 0.2; eta.ka ~ 0.1; eta.cl ~ 0.05 })
    model({ ka <- mix(exp(tka1 + eta.ka), p1, exp(tka2 + eta.ka))
            cl <- mix(exp(tcl1 + eta.cl), p1, exp(tcl2 + eta.cl))
            v <- exp(tv); cp <- linCmt(); cp ~ add(add.sd) })
  }
  cl <- .rpemClassify(rxode2::rxode2(mod))
  expect_equal(cl$mix$nParam, 2L)
  expect_equal(cl$mix$K, 2L)
  expect_equal(cl$nEta, 2L)
  expect_equal(dim(cl$mix$etaForComp), c(2L, 2L))          # nParam x K
  expect_equal(dim(cl$mix$muCompIdx), c(2L, 2L))
})

test_that("mix() calls with mismatched probabilities are rejected", {
  # rxode2 enforces one shared latent class at UI construction (before RPEM sees it);
  # .rpemMixInfo keeps a defensive check for the same invariant.
  bad <- function() {
    ini({ tka1 <- 0.2; tka2 <- 1.2; tcl1 <- 0.8; tcl2 <- 1.5; tv <- fix(3.45)
          p1 <- 0.5; q1 <- 0.5; add.sd <- 0.2; eta.ka ~ 0.1; eta.cl ~ 0.05 })
    model({ ka <- mix(exp(tka1 + eta.ka), p1, exp(tka2 + eta.ka))
            cl <- mix(exp(tcl1 + eta.cl), q1, exp(tcl2 + eta.cl))   # different probability
            v <- exp(tv); cp <- linCmt(); cp ~ add(add.sd) })
  }
  expect_error(rxode2::rxode2(bad), "probabilities in a mixture must match")
})

test_that("RPEM recovers a two-parameter mixture (ka and cl)", {
  skip_on_cran()
  skip_on_ci()  # heavy: multi-iteration mixture EM

  sim <- rxode2::rxode2({ ka <- exp(tka + eka); cl <- exp(tcl + ecl); v <- exp(tv); cp <- linCmt() })
  set.seed(42); nsub <- 150L; obsT <- seq(0.5, 24, by = 2)
  tk <- c(0.0, 1.4); tc <- c(0.9, 1.4); w1 <- 0.6
  dat <- do.call(rbind, lapply(seq_len(nsub), function(i) {
    k <- if (stats::runif(1) < w1) 1L else 2L
    ev <- rxode2::et(amt = 100, cmt = "depot"); ev <- rxode2::et(ev, obsT)
    s <- rxode2::rxSolve(sim, params = c(tka = tk[k], tcl = tc[k], tv = 3.45,
                                         eka = stats::rnorm(1, 0, sqrt(0.1)),
                                         ecl = stats::rnorm(1, 0, sqrt(0.05))),
                         events = ev, returnType = "data.frame", addDosing = FALSE)
    d <- as.data.frame(ev); d$id <- i; o <- d$evid == 0; d$DV <- 0
    d$DV[o] <- s$cp + stats::rnorm(nrow(s), 0, 0.1)
    d
  }))
  mod <- function() {
    ini({ tka1 <- 0.2; tka2 <- 1.2; tcl1 <- 0.8; tcl2 <- 1.5; tv <- fix(3.45)
          p1 <- 0.5; add.sd <- 0.2; eta.ka ~ 0.1; eta.cl ~ 0.05 })
    model({ ka <- mix(exp(tka1 + eta.ka), p1, exp(tka2 + eta.ka))
            cl <- mix(exp(tcl1 + eta.cl), p1, exp(tcl2 + eta.cl))
            v <- exp(tv); cp <- linCmt(); cp ~ add(add.sd) })
  }
  rf <- .rpemFit(rxode2::rxode2(mod), dat,
                 rpemControl(nGauss = 400L, nMH = 80000L, mhBurn = 8000L,
                             niter = 30L, collect = 10L, seed = 7L, cores = 4L))
  # muK is nParam x K: row 1 = ka (tka1, tka2), row 2 = cl (tcl1, tcl2)
  expect_equal(dim(rf$mix$muK), c(2L, 2L))
  expect_equal(unname(rf$mix$muK[1, 1]), 0.0, tolerance = 0.2)   # tka comp 1
  expect_equal(unname(rf$mix$muK[1, 2]), 1.4, tolerance = 0.2)   # tka comp 2
  expect_equal(unname(rf$mix$muK[2, 1]), 0.9, tolerance = 0.2)   # tcl comp 1
  expect_equal(unname(rf$mix$muK[2, 2]), 1.4, tolerance = 0.2)   # tcl comp 2
  expect_equal(unname(rf$mix$w[1]), 0.6, tolerance = 0.15)
  expect_equal(rf$addSd, 0.1, tolerance = 0.06)
})
