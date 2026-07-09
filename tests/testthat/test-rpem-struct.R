# Non-mu-referenced structural fixed effects: a parameter with no random effect
# (here tcl) is estimated by the numeric re-solve M-step (rpemMstepBeta), while a
# fix()ed parameter (tv) is held (design/rpem/05, D6/D20).  Must match FOCEI.

test_that("RPEM estimates a non-mu-ref structural parameter (matches FOCEI)", {
  skip_on_cran()
  skip_on_ci()  # heavy: FOCEI fit + re-solving numeric M-step

  struct <- rxode2::rxode2({ ka <- exp(tka + eta); cl <- exp(tcl); v <- exp(tv); cp <- linCmt() })
  set.seed(9); nsub <- 40L; obsT <- seq(0.5, 24, by = 1.5); trueAdd <- 0.1
  etasTrue <- rnorm(nsub, 0, sqrt(0.3))
  dat <- do.call(rbind, lapply(seq_len(nsub), function(i) {
    ev <- rxode2::et(amt = 100, cmt = "depot"); ev <- rxode2::et(ev, obsT)
    s <- rxode2::rxSolve(struct, params = c(tka = 0.45, tcl = 1, tv = 3.45, eta = etasTrue[i]),
                         events = ev, returnType = "data.frame", addDosing = FALSE)
    d <- as.data.frame(ev); d$id <- i; d$DV <- 0; o <- d$evid == 0
    d$DV[o] <- s$cp + rnorm(nrow(s), 0, trueAdd)
    d
  }))

  # tcl estimated (structural, no eta, no fix); tv fixed
  mod <- function() {
    ini({ tka <- 0.3; tcl <- 0.5; tv <- fix(3.45); add.sd <- 0.2; eta.ka ~ 1.0 })
    model({ ka <- exp(tka + eta.ka); cl <- exp(tcl); v <- exp(tv); cp <- linCmt(); cp ~ add(add.sd) })
  }
  fit <- suppressMessages(nlmixr2(mod, dat, est = "focei",
                                  control = foceiControl(print = 0, calcTables = FALSE)))
  fTka <- fit$parFixedDf["tka", "Estimate"]
  fTcl <- fit$parFixedDf["tcl", "Estimate"]; fAdd <- fit$parFixedDf["add.sd", "Estimate"]

  rmod <- function() {
    ini({ tka <- 0.3; tcl <- 0.5; tv <- fix(3.45); add.sd <- 0.2; eta.ka ~ 0.6 })
    model({ ka <- exp(tka + eta.ka); cl <- exp(tcl); v <- exp(tv); cp <- linCmt(); cp ~ add(add.sd) })
  }
  ui <- rxode2::rxode2(rmod)
  cl <- .rpemClassify(ui)
  expect_equal(cl$fixNames, "tv")
  expect_equal(cl$thetaNames[cl$structIdx + 1L], "tcl")

  rf <- .rpemFit(ui, dat, rpemControl(nGauss = 200L, nMH = 60000L, mhBurn = 6000L,
                                      niter = 25L, collect = 10L, seed = 400L))
  expect_equal(unname(rf$mu["tka"]), fTka, tolerance = 0.06)
  expect_equal(unname(rf$struct["tcl"]), fTcl, tolerance = 0.05)
  expect_equal(rf$addSd, fAdd, tolerance = 0.04)
})
