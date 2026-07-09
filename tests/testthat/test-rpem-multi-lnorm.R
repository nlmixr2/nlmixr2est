# Per-endpoint lognormal: the common PK-lognormal + PD-additive combination.  One
# endpoint uses lnorm (additive on the log scale), the other additive; each residual
# is updated over just its own observations (design/rpem/05).  Must match FOCEI.

test_that("RPEM supports lognormal on one endpoint, additive on another (matches FOCEI)", {
  skip_on_cran()
  skip_on_ci()  # heavy: FOCEI fit + multi-iteration RPEM loop

  struct <- rxode2::rxode2({ ka <- exp(tka + eta); cl <- exp(tcl); v <- exp(tv); cp <- linCmt() })
  set.seed(7); nsub <- 40L; obsT <- seq(0.5, 24, by = 2); te0 <- 0.7; lnTrue <- 0.15; effTrue <- 0.5
  etasTrue <- rnorm(nsub, 0, sqrt(0.3))
  dat <- do.call(rbind, lapply(seq_len(nsub), function(i) {
    ev <- rxode2::et(amt = 100, cmt = "depot"); ev <- rxode2::et(ev, obsT)
    s <- rxode2::rxSolve(struct, params = c(tka = 0.45, tcl = 1, tv = 3.45, eta = etasTrue[i]),
                         events = ev, returnType = "data.frame", addDosing = FALSE)
    cp <- s$cp; eff <- cp * exp(te0)
    dose <- data.frame(id = i, time = 0, dvid = 1L, DV = 0, evid = 1L, amt = 100, cmt = 1L)
    d1 <- data.frame(id = i, time = obsT, dvid = 1L, DV = cp * exp(rnorm(length(cp), 0, lnTrue)),
                     evid = 0L, amt = 0, cmt = 1L)
    d2 <- data.frame(id = i, time = obsT, dvid = 2L, DV = eff + rnorm(length(eff), 0, effTrue),
                     evid = 0L, amt = 0, cmt = 2L)
    rbind(dose, d1, d2)
  }))
  dat <- dat[order(dat$id, dat$time, dat$dvid), ]

  mod <- function() {
    ini({ tka <- 0.3; tcl <- fix(1.0); tv <- fix(3.45); te0 <- fix(0.7)
          prop.sd <- 0.2; eff.sd <- 0.3; eta.ka ~ 1.0 })
    model({ ka <- exp(tka + eta.ka); cl <- exp(tcl); v <- exp(tv); cp <- linCmt(); eff <- cp * exp(te0)
            cp ~ lnorm(prop.sd); eff ~ add(eff.sd) })
  }
  rmod <- function() {
    ini({ tka <- 0.3; tcl <- fix(1.0); tv <- fix(3.45); te0 <- fix(0.7)
          prop.sd <- 0.2; eff.sd <- 0.3; eta.ka ~ 0.6 })
    model({ ka <- exp(tka + eta.ka); cl <- exp(tcl); v <- exp(tv); cp <- linCmt(); eff <- cp * exp(te0)
            cp ~ lnorm(prop.sd); eff ~ add(eff.sd) })
  }
  fit <- suppressMessages(nlmixr2(mod, dat, est = "focei",
                                  control = foceiControl(print = 0, calcTables = FALSE)))
  fTka <- fit$parFixedDf["tka", "Estimate"]
  fLn <- fit$parFixedDf["prop.sd", "Estimate"]; fEff <- fit$parFixedDf["eff.sd", "Estimate"]

  ui <- rxode2::rxode2(rmod)
  cl <- .rpemClassify(ui)
  expect_equal(cl$endpt$errType, c(6L, 0L))   # cp lognormal, eff additive
  rf <- .rpemFit(ui, dat, rpemControl(nGauss = 300L, nMH = 80000L, mhBurn = 8000L,
                                      niter = 30L, collect = 12L, seed = 300L, cores = 4L))
  expect_equal(unname(rf$mu["tka"]), fTka, tolerance = 0.08)
  expect_equal(rf$endptSd[1], fLn, tolerance = 0.03)    # cp lnorm sd
  expect_equal(rf$endptSd[2], fEff, tolerance = 0.05)   # eff add.sd
})
