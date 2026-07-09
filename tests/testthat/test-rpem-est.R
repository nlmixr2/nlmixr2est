# est="rpem" dispatch: nlmixr2(model, data, est="rpem") must resolve the control,
# run the E-M loop, and return a fit object with sensible estimates (K=1 core).

test_that("est='rpem' dispatch runs end-to-end (K=1)", {
  skip_on_cran()
  skip_on_ci()  # heavy: multi-iteration RPEM loop

  struct <- rxode2::rxode2({ ka <- exp(tka+eta); cl <- exp(tcl); v <- exp(tv); cp <- linCmt() })
  set.seed(3)
  nsub <- 25L; addSd <- 0.15; obsT <- seq(0.25, 24, by=1.5)
  etasTrue <- rnorm(nsub, 0, sqrt(0.3))
  dat <- do.call(rbind, lapply(seq_len(nsub), function(i) {
    ev <- rxode2::et(amt=100, cmt="depot"); ev <- rxode2::et(ev, obsT)
    s <- rxode2::rxSolve(struct, params=c(tka=0.45, tcl=1, tv=3.45, eta=etasTrue[i]),
                         events=ev, returnType="data.frame", addDosing=FALSE)
    d <- as.data.frame(ev); d$id <- i; d$DV <- 0; o <- d$evid == 0
    d$DV[o] <- s$cp + rnorm(nrow(s), 0, addSd); d
  }))

  rmod <- function() {
    ini({ tka <- 0.3; tcl <- 1.0; tv <- 3.45; add.sd <- 0.3; eta.ka ~ 0.6 })
    model({ ka <- exp(tka + eta.ka); cl <- exp(tcl); v <- exp(tv); cp <- linCmt(); cp ~ add(add.sd) })
  }

  fit <- nlmixr2(rmod, dat, est="rpem",
                 control=rpemControl(nGauss=200L, nMH=50000L, mhBurn=5000L,
                                     niter=20L, collect=8L, seed=100L))

  # Plumbing smoke test: dispatch runs and returns finite, plausible estimates.
  # Precise FOCEI-recovery is covered by test-rpem-fit.R.
  expect_s3_class(fit, "nlmixr2rpem")
  expect_true(is.finite(unname(fit$mu["tka"])))
  expect_gt(unname(fit$mu["tka"]), 0.1);  expect_lt(unname(fit$mu["tka"]), 1.2)
  expect_gt(unname(fit$omega["eta.ka"]), 0.05); expect_lt(unname(fit$omega["eta.ka"]), 0.8)
  expect_gt(fit$addSd, 0.05); expect_lt(fit$addSd, 0.3)
})
