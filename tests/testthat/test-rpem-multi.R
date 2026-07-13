# Multiple endpoints (mirror SAEM): the E-step computes the joint multi-endpoint
# likelihood and the M-step updates each endpoint's residual separately over its
# own observations (design/rpem/05).  Must match FOCEI.

test_that("RPEM supports multiple endpoints with per-endpoint residuals (matches FOCEI)", {
  skip_on_cran()

  struct <- rxode2::rxode2({ ka <- exp(tka + eta); cl <- exp(tcl); v <- exp(tv); cp <- linCmt() })
  set.seed(5); nsub <- 40L; obsT <- seq(0.5, 24, by = 2); trueAdd <- 0.1; trueEff <- 0.5; te0 <- 0.7
  etasTrue <- rnorm(nsub, 0, sqrt(0.3))
  dat <- do.call(rbind, lapply(seq_len(nsub), function(i) {
    ev <- rxode2::et(amt = 100, cmt = "depot"); ev <- rxode2::et(ev, obsT)
    s <- rxode2::rxSolve(struct, params = c(tka = 0.45, tcl = 1, tv = 3.45, eta = etasTrue[i]),
                         events = ev, returnType = "data.frame", addDosing = FALSE)
    cp <- s$cp; eff <- cp * exp(te0)
    dose <- data.frame(id = i, time = 0, dvid = 1L, DV = 0, evid = 1L, amt = 100, cmt = 1L)
    d1 <- data.frame(id = i, time = obsT, dvid = 1L, DV = cp + rnorm(length(cp), 0, trueAdd),
                     evid = 0L, amt = 0, cmt = 1L)
    d2 <- data.frame(id = i, time = obsT, dvid = 2L, DV = eff + rnorm(length(eff), 0, trueEff),
                     evid = 0L, amt = 0, cmt = 2L)
    rbind(dose, d1, d2)
  }))
  dat <- dat[order(dat$id, dat$time, dat$dvid), ]

  mod <- function() {
    ini({ tka <- 0.3; tcl <- fix(1.0); tv <- fix(3.45); te0 <- fix(0.7); add.sd <- 0.3; eff.sd <- 0.3; eta.ka ~ 1.0 })
    model({ ka <- exp(tka + eta.ka); cl <- exp(tcl); v <- exp(tv); cp <- linCmt(); eff <- cp * exp(te0)
            cp ~ add(add.sd); eff ~ add(eff.sd) })
  }
  fit <- suppressMessages(nlmixr2(mod, dat, est = "focei",
                                  control = foceiControl(print = 0, calcTables = FALSE)))
  fTka <- fit$parFixedDf["tka", "Estimate"]
  fAdd <- fit$parFixedDf["add.sd", "Estimate"]; fEff <- fit$parFixedDf["eff.sd", "Estimate"]

  rmod <- function() {
    ini({ tka <- 0.3; tcl <- fix(1.0); tv <- fix(3.45); te0 <- fix(0.7); add.sd <- 0.3; eff.sd <- 0.3; eta.ka ~ 0.6 })
    model({ ka <- exp(tka + eta.ka); cl <- exp(tcl); v <- exp(tv); cp <- linCmt(); eff <- cp * exp(te0)
            cp ~ add(add.sd); eff ~ add(eff.sd) })
  }
  ui <- rxode2::rxode2(rmod)
  rf <- .rpemFit(ui, dat, rpemControl(nGauss = 300L, nMH = 80000L, mhBurn = 8000L,
                                      niter = 30L, collect = 12L, seed = 300L))

  cl <- .rpemClassify(ui)
  expect_equal(cl$errName, "multiEndpoint")
  expect_equal(cl$endpt$nEndpt, 2L)
  expect_equal(unname(rf$mu["tka"]), fTka, tolerance = 0.08)
  expect_equal(rf$endptSd[1], fAdd, tolerance = 0.05)  # cp additive
  expect_equal(rf$endptSd[2], fEff, tolerance = 0.05)  # eff additive
})

test_that("RPEM multi-endpoint per-endpoint residual M-step recovers each endpoint's error", {
  skip_on_cran()

  struct <- rxode2::rxode2({ ka <- exp(tka + eta); cl <- exp(tcl); v <- exp(tv); cp <- linCmt() })
  set.seed(7); nsub <- 45L; obsT <- seq(1, 24, by = 3)
  etasTrue <- rnorm(nsub, 0, sqrt(0.3))
  dat <- do.call(rbind, lapply(seq_len(nsub), function(i) {
    ev <- rxode2::et(amt = 100, cmt = "depot"); ev <- rxode2::et(ev, obsT)
    s <- rxode2::rxSolve(struct, params = c(tka = 0.45, tcl = 1, tv = 3.45, eta = etasTrue[i]),
                         events = ev, returnType = "data.frame", addDosing = FALSE)
    cp <- s$cp; eff <- cp * exp(0.7)
    dose <- data.frame(id = i, time = 0, dvid = 1L, DV = 0, evid = 1L, amt = 100, cmt = 1L)
    # endpoint 1: combined (add + prop); endpoint 2: proportional
    d1 <- data.frame(id = i, time = obsT, dvid = 1L, evid = 0L, amt = 0, cmt = 1L,
                     DV = cp + rnorm(length(cp), 0, sqrt(0.05^2 + (0.1 * cp)^2)))
    d2 <- data.frame(id = i, time = obsT, dvid = 2L, evid = 0L, amt = 0, cmt = 2L,
                     DV = eff * (1 + rnorm(length(eff), 0, 0.12)))
    rbind(dose, d1, d2)
  }))
  dat <- dat[order(dat$id, dat$time, dat$dvid), ]
  rmod <- function() {
    ini({ tka <- 0.3; tcl <- fix(1.0); tv <- fix(3.45); te0 <- fix(0.7)
          add.sd <- 0.1; prop.sd <- 0.15; eff.sd <- 0.2; eta.ka ~ 0.6 })
    model({ ka <- exp(tka + eta.ka); cl <- exp(tcl); v <- exp(tv); cp <- linCmt(); eff <- cp * exp(te0)
            cp ~ add(add.sd) + prop(prop.sd); eff ~ prop(eff.sd) })
  }
  ui <- rxode2::rxUiDecompress(rxode2::rxode2(rmod))
  expect_equal(.rpemClassify(ui)$errType, 5L)
  rfC <- .rpemFit(ui, dat, rpemControl(nGauss = 300L, nMH = 60000L, mhBurn = 6000L, niter = 25L,
                                       collect = 10L, seed = 42L, cores = 4L))
  # the per-endpoint residual M-step (combined ep1 + proportional ep2) recovers tka and each
  # endpoint's residual in the right ballpark (wide bands: the add / prop split of a combined
  # error is weakly identified)
  expect_equal(unname(rfC$mu["tka"]), 0.45, tolerance = 0.12)
  expect_gt(rfC$endptSd[1], 0.01); expect_lt(rfC$endptSd[1], 0.09)     # ep1 add ~0.05
  expect_equal(rfC$endptProp[1], 0.10, tolerance = 0.3)               # ep1 prop ~0.10
  expect_equal(rfC$endptSd[2], 0.12, tolerance = 0.25)               # ep2 prop ~0.12
})
