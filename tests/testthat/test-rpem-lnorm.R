# Lognormal residual error (lnorm): additive error on the log scale (a fixed log
# transform whose Jacobian is constant in the parameters, so it cancels in the
# E-step / MH).  The residual sd is the RMSE of the log residuals (design/rpem/05).

test_that("RPEM supports single-endpoint lognormal error (matches FOCEI)", {
  skip_on_cran()
  skip_on_ci()  # heavy: FOCEI fit + multi-iteration RPEM loop

  struct <- rxode2::rxode2({ ka <- exp(tka + eta); cl <- exp(tcl); v <- exp(tv); cp <- linCmt() })
  set.seed(7); nsub <- 40L; obsT <- seq(0.5, 24, by = 1.5); lnTrue <- 0.15
  etasTrue <- rnorm(nsub, 0, sqrt(0.3))
  dat <- do.call(rbind, lapply(seq_len(nsub), function(i) {
    ev <- rxode2::et(amt = 100, cmt = "depot"); ev <- rxode2::et(ev, obsT)
    s <- rxode2::rxSolve(struct, params = c(tka = 0.45, tcl = 1, tv = 3.45, eta = etasTrue[i]),
                         events = ev, returnType = "data.frame", addDosing = FALSE)
    d <- as.data.frame(ev); d$id <- i; d$DV <- 0; o <- d$evid == 0
    d$DV[o] <- s$cp * exp(rnorm(nrow(s), 0, lnTrue))
    d
  }))

  mod <- function() {
    ini({ tka <- 0.3; tcl <- fix(1.0); tv <- fix(3.45); prop.sd <- 0.2; eta.ka ~ 1.0 })
    model({ ka <- exp(tka + eta.ka); cl <- exp(tcl); v <- exp(tv); cp <- linCmt(); cp ~ lnorm(prop.sd) })
  }
  fit <- suppressMessages(nlmixr2(mod, dat, est = "focei",
                                  control = foceiControl(print = 0, calcTables = FALSE)))
  fTka <- fit$parFixedDf["tka", "Estimate"]; fSd <- fit$parFixedDf["prop.sd", "Estimate"]

  rmod <- function() {
    ini({ tka <- 0.3; tcl <- fix(1.0); tv <- fix(3.45); prop.sd <- 0.2; eta.ka ~ 0.6 })
    model({ ka <- exp(tka + eta.ka); cl <- exp(tcl); v <- exp(tv); cp <- linCmt(); cp ~ lnorm(prop.sd) })
  }
  ui <- rxode2::rxode2(rmod)
  expect_equal(.rpemClassify(ui)$errName, "lnorm")
  rf <- .rpemFit(ui, dat, rpemControl(nGauss = 300L, nMH = 80000L, mhBurn = 8000L,
                                      niter = 30L, collect = 12L, seed = 100L, cores = 4L))
  expect_equal(unname(rf$mu["tka"]), fTka, tolerance = 0.06)
  expect_equal(rf$addSd, fSd, tolerance = 0.03)   # addSd slot holds the lnorm sd
})
