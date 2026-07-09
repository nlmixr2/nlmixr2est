# Combined residual error (add + prop): the combined variance add^2 + prop^2 cp^2
# has no closed form, so RPEM's M-step Newton-optimizes (add.sd, prop.sd) over the
# accepted samples' per-obs cp/residual and must match FOCEI (design/rpem/05, D20).

test_that("RPEM supports combined (add + prop) error (matches FOCEI)", {
  skip_on_cran()
  skip_on_ci()  # heavy: FOCEI fit + multi-iteration RPEM loop

  struct <- rxode2::rxode2({ ka <- exp(tka + eta); cl <- exp(tcl); v <- exp(tv); cp <- linCmt() })
  set.seed(42)
  nsub <- 40L; trueTka <- 0.45; addSd <- 0.2; propSd <- 0.1; obsT <- seq(0.25, 24, by = 1.5)
  etasTrue <- rnorm(nsub, 0, sqrt(0.3))
  dat <- do.call(rbind, lapply(seq_len(nsub), function(i) {
    ev <- rxode2::et(amt = 100, cmt = "depot"); ev <- rxode2::et(ev, obsT)
    s <- rxode2::rxSolve(struct, params = c(tka = trueTka, tcl = 1, tv = 3.45, eta = etasTrue[i]),
                         events = ev, returnType = "data.frame", addDosing = FALSE,
                         atol = 1e-8, rtol = 1e-8)
    d <- as.data.frame(ev); d$id <- i; d$DV <- 0; o <- d$evid == 0
    d$DV[o] <- s$cp + rnorm(nrow(s), 0, sqrt(addSd^2 + (propSd * s$cp)^2))  # combined noise
    d
  }))

  mod <- function() {
    ini({ tka <- 0.3; tcl <- fix(1.0); tv <- fix(3.45); add.sd <- 0.3; prop.sd <- 0.2; eta.ka ~ 1.0 })
    model({ ka <- exp(tka + eta.ka); cl <- exp(tcl); v <- exp(tv); cp <- linCmt(); cp ~ add(add.sd) + prop(prop.sd) })
  }
  fit <- suppressMessages(nlmixr2(mod, dat, est = "focei",
                                  control = foceiControl(print = 0, calcTables = FALSE)))
  fTka <- fit$parFixedDf["tka", "Estimate"]
  fAdd <- fit$parFixedDf["add.sd", "Estimate"]; fProp <- fit$parFixedDf["prop.sd", "Estimate"]

  rmod <- function() {
    ini({ tka <- 0.3; tcl <- 1.0; tv <- 3.45; add.sd <- 0.3; prop.sd <- 0.2; eta.ka ~ 0.6 })
    model({ ka <- exp(tka + eta.ka); cl <- exp(tcl); v <- exp(tv); cp <- linCmt(); cp ~ add(add.sd) + prop(prop.sd) })
  }
  ui <- rxode2::rxode2(rmod)
  rf <- .rpemFit(ui, dat, rpemControl(nGauss = 300L, nMH = 80000L, mhBurn = 8000L,
                                      niter = 30L, collect = 12L, seed = 3000L))

  expect_equal(.rpemClassify(ui)$errName, "add+prop")
  expect_equal(unname(rf$mu["tka"]), fTka, tolerance = 0.05)
  expect_equal(rf$addSd, fAdd, tolerance = 0.05)
  expect_equal(rf$propSd, fProp, tolerance = 0.03)
})
