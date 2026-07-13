# Iteration printing + parameter history for est="rpem": RPEM uses the shared scale.h
# iteration-print machinery (like saem / focei / vae) to print the population estimates each
# iteration (with a back-transformed row and the algorithm-phase label) and to capture the
# parameter-history walk, which is installed as standard parHistData on the fit object.

.phData <- function(seed = 8L, nsub = 25L) {
  struct <- rxode2::rxode2({ ka <- exp(tka + eta); cl <- exp(tcl); v <- exp(tv); cp <- linCmt() })
  set.seed(seed); etas <- rnorm(nsub, 0, sqrt(0.3)); obsT <- seq(0.5, 24, by = 3)
  do.call(rbind, lapply(seq_len(nsub), function(i) {
    ev <- rxode2::et(amt = 100, cmt = "depot"); ev <- rxode2::et(ev, obsT)
    s <- rxode2::rxSolve(struct, params = c(tka = 0.45, tcl = 1, tv = 3.45, eta = etas[i]),
                         events = ev, returnType = "data.frame", addDosing = FALSE)
    d <- as.data.frame(ev); d$id <- i; o <- d$evid == 0; d$DV <- 0
    d$DV[o] <- s$cp + rnorm(nrow(s), 0, 0.1); d
  }))
}

.phMod <- function() {
  ini({ tka <- 0.3; tcl <- fix(1.0); tv <- fix(3.45); add.sd <- 0.2; eta.ka ~ 0.6 })
  model({ ka <- exp(tka + eta.ka); cl <- exp(tcl); v <- exp(tv); cp <- linCmt(); cp ~ add(add.sd) })
}

test_that("rpemControl builds an iterPrintControl and defaults print to 1 (live trace)", {
  ctl <- rpemControl()
  expect_equal(ctl$print, 1L)                    # saem/focei/vae-style live iteration trace
  expect_s3_class(ctl$iterPrintControl, "iterPrintControl")
  expect_equal(rpemControl(print = 0L)$print, 0L)  # opt-in silent capture still available
})

test_that("RPEM captures the parameter-history walk (one row per iteration + a phase)", {
  skip_on_cran()

  ui <- rxode2::rxUiDecompress(rxode2::rxode2(.phMod))
  rf <- .rpemFit(ui, .phData(), rpemControl(nGauss = 200L, nMH = 30000L, mhBurn = 3000L,
                                            niter = 12L, collect = 4L, seed = 1L, cores = 4L))
  expect_s3_class(rf$parHist, "data.frame")
  # standard parHistData columns: iter, type, objf, then the estimated params + omega
  expect_true(all(c("iter", "type", "objf") %in% names(rf$parHist)))
  expect_true("tka" %in% names(rf$parHist))
  expect_true("o(eta.ka)" %in% names(rf$parHist))
  expect_equal(max(rf$parHist$iter), 12L)                 # one iter index per EM iteration
})

test_that("a full est=rpem fit exposes fit$parHist and fit$parHistStacked", {
  skip_on_cran()

  f <- suppressMessages(suppressWarnings(nlmixr2est::nlmixr2(.phMod, .phData(), est = "rpem",
    control = rpemControl(nGauss = 200L, nMH = 30000L, mhBurn = 3000L, niter = 12L,
                          collect = 4L, seed = 1L, cores = 4L))))
  ph <- f$parHist
  expect_s3_class(ph, "data.frame")
  expect_equal(nrow(ph), 12L)
  expect_true("tka" %in% names(ph))
  phs <- f$parHistStacked
  expect_s3_class(phs, "data.frame")
  expect_gt(nrow(phs), 0)
})

test_that("RPEM iteration printing (print > 0) runs and captures the walk", {
  skip_on_cran()

  ui <- rxode2::rxUiDecompress(rxode2::rxode2(.phMod))
  # print = 2 displays the iteration table (C-level RSprintf) every 2 iterations; the walk
  # is captured regardless.  (The rendered table is verified by eye; here we assert the
  # display path runs without error and the parameter history is produced.)
  rf <- suppressMessages(
    .rpemFit(ui, .phData(), rpemControl(nGauss = 150L, nMH = 20000L, mhBurn = 2000L,
                                        niter = 10L, collect = 3L, seed = 1L, cores = 4L,
                                        print = 2L)))
  expect_s3_class(rf$parHist, "data.frame")
  expect_equal(max(rf$parHist$iter), 10L)
})

test_that("RPEM mixture fits capture the per-component parameter history", {
  skip_on_cran()

  sim <- rxode2::rxode2({ ka <- exp(tka + eka); cl <- exp(tcl); v <- exp(tv); cp <- linCmt() })
  set.seed(52); nsub <- 100L; obsT <- seq(0.5, 24, by = 2)
  dat <- do.call(rbind, lapply(seq_len(nsub), function(i) {
    k <- if (stats::runif(1) < 0.5) 1L else 2L
    ev <- rxode2::et(amt = 100, cmt = "depot"); ev <- rxode2::et(ev, obsT)
    s <- rxode2::rxSolve(sim, params = c(tka = c(0.1, 1.3)[k], tcl = 1, tv = 3.45,
                                         eka = stats::rnorm(1, 0, sqrt(0.3))),
                         events = ev, returnType = "data.frame", addDosing = FALSE)
    d <- as.data.frame(ev); d$id <- i; o <- d$evid == 0; d$DV <- 0
    d$DV[o] <- s$cp + rnorm(sum(o), 0, 0.1); d
  }))
  mod <- function() {
    ini({ tka1 <- 0.0; tka2 <- 1.1; tcl <- fix(1.0); tv <- fix(3.45); p1 <- 0.5
          add.sd <- 0.1; eta.ka ~ 0.3 })
    model({ ka <- mix(exp(tka1 + eta.ka), p1, exp(tka2 + eta.ka)); cl <- exp(tcl)
            v <- exp(tv); cp <- linCmt(); cp ~ add(add.sd) })
  }
  ui <- rxode2::rxUiDecompress(rxode2::rxode2(mod))
  rf <- .rpemFit(ui, dat, rpemControl(nGauss = 200L, nMH = 40000L, mhBurn = 4000L,
                                      niter = 12L, collect = 4L, seed = 7L, cores = 4L))
  expect_s3_class(rf$parHist, "data.frame")
  # per-component typical values + the mixture weight are in the walk
  expect_true(all(c("tka1", "tka2", "p1", "o(eta.ka)") %in% names(rf$parHist)))
  expect_equal(max(rf$parHist$iter), 12L)
})
