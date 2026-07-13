# Power residual error: variance (prop.sd * cp^power)^2 with an estimated exponent.
# RPEM profiles out the scale and golden-section optimizes the exponent (design/
# rpem/05, D20).  Must match FOCEI.

test_that("RPEM supports power error with an estimated exponent (matches FOCEI)", {
  skip_on_cran()

  struct <- rxode2::rxode2({ ka <- exp(tka + eta); cl <- exp(tcl); v <- exp(tv); cp <- linCmt() })
  set.seed(11); nsub <- 40L; obsT <- seq(0.5, 24, by = 1.5); truePropSd <- 0.15; truePow <- 0.75
  etasTrue <- rnorm(nsub, 0, sqrt(0.3))
  dat <- do.call(rbind, lapply(seq_len(nsub), function(i) {
    ev <- rxode2::et(amt = 100, cmt = "depot"); ev <- rxode2::et(ev, obsT)
    s <- rxode2::rxSolve(struct, params = c(tka = 0.45, tcl = 1, tv = 3.45, eta = etasTrue[i]),
                         events = ev, returnType = "data.frame", addDosing = FALSE)
    d <- as.data.frame(ev); d$id <- i; d$DV <- 0; o <- d$evid == 0
    d$DV[o] <- s$cp + rnorm(nrow(s), 0, truePropSd * (abs(s$cp)^truePow))
    d
  }))

  mod <- function() {
    ini({ tka <- 0.3; tcl <- fix(1.0); tv <- fix(3.45); prop.sd <- 0.2; pw <- 0.5; eta.ka ~ 1.0 })
    model({ ka <- exp(tka + eta.ka); cl <- exp(tcl); v <- exp(tv); cp <- linCmt(); cp ~ pow(prop.sd, pw) })
  }
  fit <- suppressMessages(nlmixr2(mod, dat, est = "focei",
                                  control = foceiControl(print = 0, calcTables = FALSE)))
  fTka <- fit$parFixedDf["tka", "Estimate"]
  fProp <- fit$parFixedDf["prop.sd", "Estimate"]; fPw <- fit$parFixedDf["pw", "Estimate"]

  rmod <- function() {
    ini({ tka <- 0.3; tcl <- fix(1.0); tv <- fix(3.45); prop.sd <- 0.2; pw <- 0.5; eta.ka ~ 0.6 })
    model({ ka <- exp(tka + eta.ka); cl <- exp(tcl); v <- exp(tv); cp <- linCmt(); cp ~ pow(prop.sd, pw) })
  }
  ui <- rxode2::rxode2(rmod)
  # power-model location needs more MC precision (late-time low-cp points dominate)
  rf <- .rpemFit(ui, dat, rpemControl(nGauss = 400L, nMH = 100000L, mhBurn = 10000L,
                                      niter = 40L, collect = 15L, seed = 200L))

  expect_equal(.rpemClassify(ui)$errName, "pow")
  # tolerance is relative (all.equal); tka is small (~0.23) so allow ~15%
  expect_equal(unname(rf$mu["tka"]), fTka, tolerance = 0.15)
  expect_equal(rf$addSd, fProp, tolerance = 0.05)   # addSd slot holds the power scale (prop.sd)
  expect_equal(rf$power, fPw, tolerance = 0.1)
})
