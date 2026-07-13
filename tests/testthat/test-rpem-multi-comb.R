# Per-endpoint combined error: one endpoint uses combined (add + prop) error, the
# other additive.  The per-endpoint M-step runs the guarded combined optimizer over
# just that endpoint's observations (design/rpem/05).  Must match FOCEI.

test_that("RPEM supports combined error on one of several endpoints (matches FOCEI)", {
  skip_on_cran()

  struct <- rxode2::rxode2({ ka <- exp(tka + eta); cl <- exp(tcl); v <- exp(tv); cp <- linCmt() })
  set.seed(5); nsub <- 40L; obsT <- seq(0.5, 24, by = 2); te0 <- 0.7
  addTrue <- 0.1; propTrue <- 0.08; effTrue <- 0.5
  etasTrue <- rnorm(nsub, 0, sqrt(0.3))
  dat <- do.call(rbind, lapply(seq_len(nsub), function(i) {
    ev <- rxode2::et(amt = 100, cmt = "depot"); ev <- rxode2::et(ev, obsT)
    s <- rxode2::rxSolve(struct, params = c(tka = 0.45, tcl = 1, tv = 3.45, eta = etasTrue[i]),
                         events = ev, returnType = "data.frame", addDosing = FALSE)
    cp <- s$cp; eff <- cp * exp(te0)
    dose <- data.frame(id = i, time = 0, dvid = 1L, DV = 0, evid = 1L, amt = 100, cmt = 1L)
    d1 <- data.frame(id = i, time = obsT, dvid = 1L,
                     DV = cp + rnorm(length(cp), 0, sqrt(addTrue^2 + (propTrue * cp)^2)),
                     evid = 0L, amt = 0, cmt = 1L)
    d2 <- data.frame(id = i, time = obsT, dvid = 2L, DV = eff + rnorm(length(eff), 0, effTrue),
                     evid = 0L, amt = 0, cmt = 2L)
    rbind(dose, d1, d2)
  }))
  dat <- dat[order(dat$id, dat$time, dat$dvid), ]

  mod <- function() {
    ini({ tka <- 0.3; tcl <- fix(1.0); tv <- fix(3.45); te0 <- fix(0.7)
          add.sd <- 0.3; prop.sd <- 0.2; eff.sd <- 0.3; eta.ka ~ 1.0 })
    model({ ka <- exp(tka + eta.ka); cl <- exp(tcl); v <- exp(tv); cp <- linCmt(); eff <- cp * exp(te0)
            cp ~ add(add.sd) + prop(prop.sd); eff ~ add(eff.sd) })
  }
  rmod <- function() {
    ini({ tka <- 0.3; tcl <- fix(1.0); tv <- fix(3.45); te0 <- fix(0.7)
          add.sd <- 0.3; prop.sd <- 0.2; eff.sd <- 0.3; eta.ka ~ 0.6 })
    model({ ka <- exp(tka + eta.ka); cl <- exp(tcl); v <- exp(tv); cp <- linCmt(); eff <- cp * exp(te0)
            cp ~ add(add.sd) + prop(prop.sd); eff ~ add(eff.sd) })
  }
  fit <- suppressMessages(nlmixr2(mod, dat, est = "focei",
                                  control = foceiControl(print = 0, calcTables = FALSE)))
  fTka <- fit$parFixedDf["tka", "Estimate"]
  fAdd <- fit$parFixedDf["add.sd", "Estimate"]; fProp <- fit$parFixedDf["prop.sd", "Estimate"]
  fEff <- fit$parFixedDf["eff.sd", "Estimate"]

  ui <- rxode2::rxode2(rmod)
  cl <- .rpemClassify(ui)
  expect_equal(cl$endpt$errType, c(2L, 0L))   # cp combined, eff additive
  rf <- .rpemFit(ui, dat, rpemControl(nGauss = 200L, nMH = 60000L, mhBurn = 6000L,
                                      niter = 25L, collect = 10L, seed = 300L, cores = 4L))
  expect_equal(unname(rf$mu["tka"]), fTka, tolerance = 0.1)
  expect_equal(rf$endptSd[1], fAdd, tolerance = 0.05)    # cp add.sd
  expect_equal(rf$endptProp[1], fProp, tolerance = 0.04) # cp prop.sd
  expect_equal(rf$endptSd[2], fEff, tolerance = 0.05)    # eff add.sd
})
