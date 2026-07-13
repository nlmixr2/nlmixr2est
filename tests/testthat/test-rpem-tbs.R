# Transform-both-sides (Box-Cox) with a dynamic lambda: RPEM estimates lambda +
# add.sd (transformed scale) by a golden-section profile M-step with the Jacobian
# (design/rpem/05, D20).  Must match FOCEI.

test_that("RPEM supports TBS Box-Cox with a dynamic lambda (matches FOCEI)", {
  skip_on_cran()

  struct <- rxode2::rxode2({ ka <- exp(tka + eta); cl <- exp(tcl); v <- exp(tv); cp <- linCmt() })
  set.seed(7); nsub <- 40L; obsT <- seq(0.5, 24, by = 1.5); trueLam <- 0.5; trueAdd <- 0.1
  etasTrue <- rnorm(nsub, 0, sqrt(0.3))
  bcI <- function(y, l) if (l == 0) exp(y) else (l * y + 1)^(1 / l)
  dat <- do.call(rbind, lapply(seq_len(nsub), function(i) {
    ev <- rxode2::et(amt = 100, cmt = "depot"); ev <- rxode2::et(ev, obsT)
    s <- rxode2::rxSolve(struct, params = c(tka = 0.45, tcl = 1, tv = 3.45, eta = etasTrue[i]),
                         events = ev, returnType = "data.frame", addDosing = FALSE)
    d <- as.data.frame(ev); d$id <- i; d$DV <- 0; o <- d$evid == 0
    tcp <- ((s$cp)^trueLam - 1) / trueLam
    d$DV[o] <- pmax(bcI(tcp + rnorm(nrow(s), 0, trueAdd), trueLam), 1e-3)
    d
  }))

  mod <- function() {
    ini({ tka <- 0.3; tcl <- fix(1.0); tv <- fix(3.45); add.sd <- 0.2; lambda <- 1; eta.ka ~ 1.0 })
    model({ ka <- exp(tka + eta.ka); cl <- exp(tcl); v <- exp(tv); cp <- linCmt(); cp ~ add(add.sd) + boxCox(lambda) })
  }
  fit <- suppressMessages(nlmixr2(mod, dat, est = "focei",
                                  control = foceiControl(print = 0, calcTables = FALSE)))
  fTka <- fit$parFixedDf["tka", "Estimate"]
  fAdd <- fit$parFixedDf["add.sd", "Estimate"]; fLam <- fit$parFixedDf["lambda", "Estimate"]

  rmod <- function() {
    ini({ tka <- 0.3; tcl <- fix(1.0); tv <- fix(3.45); add.sd <- 0.2; lambda <- 1; eta.ka ~ 0.6 })
    model({ ka <- exp(tka + eta.ka); cl <- exp(tcl); v <- exp(tv); cp <- linCmt(); cp ~ add(add.sd) + boxCox(lambda) })
  }
  ui <- rxode2::rxode2(rmod)
  rf <- .rpemFit(ui, dat, rpemControl(nGauss = 300L, nMH = 80000L, mhBurn = 8000L,
                                      niter = 30L, collect = 12L, seed = 100L))

  expect_equal(.rpemClassify(ui)$errName, "add+boxCox")
  expect_equal(unname(rf$mu["tka"]), fTka, tolerance = 0.06)
  expect_equal(rf$addSd, fAdd, tolerance = 0.03)
  expect_equal(rf$lambda, fLam, tolerance = 0.1)
})
