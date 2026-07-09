# Per-endpoint TBS: one endpoint uses transform-both-sides (Box-Cox with a dynamic
# lambda), the other additive.  The per-endpoint M-step profiles lambda over just
# that endpoint's observations, using the per-obs transform code cached in the
# E-step (design/rpem/05).  Validated against the simulation truth -- FOCEI aborts
# on this particular multi-endpoint Box-Cox model (a FOCEI numerical limitation
# near cp ~ 0), so it cannot serve as the reference here.

test_that("RPEM estimates TBS (Box-Cox) on one endpoint, additive on another", {
  skip_on_cran()
  skip_on_ci()  # heavy: multi-iteration RPEM loop

  struct <- rxode2::rxode2({ ka <- exp(tka + eta); cl <- exp(tcl); v <- exp(tv); cp <- linCmt() })
  set.seed(7); nsub <- 40L; obsT <- seq(0.5, 24, by = 2); te0 <- 0.7
  lamT <- 0.5; addT <- 0.1; effT <- 0.5
  etasTrue <- rnorm(nsub, 0, sqrt(0.3)); bcI <- function(y, l) if (l == 0) exp(y) else (l * y + 1)^(1 / l)
  dat <- do.call(rbind, lapply(seq_len(nsub), function(i) {
    ev <- rxode2::et(amt = 100, cmt = "depot"); ev <- rxode2::et(ev, obsT)
    s <- rxode2::rxSolve(struct, params = c(tka = 0.45, tcl = 1, tv = 3.45, eta = etasTrue[i]),
                         events = ev, returnType = "data.frame", addDosing = FALSE)
    cp <- s$cp; eff <- cp * exp(te0)
    dose <- data.frame(id = i, time = 0, dvid = 1L, DV = 0, evid = 1L, amt = 100, cmt = 1L)
    tcp <- ((cp)^lamT - 1) / lamT
    d1 <- data.frame(id = i, time = obsT, dvid = 1L,
                     DV = pmax(bcI(tcp + rnorm(length(cp), 0, addT), lamT), 1e-3),
                     evid = 0L, amt = 0, cmt = 1L)
    d2 <- data.frame(id = i, time = obsT, dvid = 2L, DV = eff + rnorm(length(eff), 0, effT),
                     evid = 0L, amt = 0, cmt = 2L)
    rbind(dose, d1, d2)
  }))
  dat <- dat[order(dat$id, dat$time, dat$dvid), ]

  rmod <- function() {
    ini({ tka <- 0.3; tcl <- fix(1.0); tv <- fix(3.45); te0 <- fix(0.7)
          add.sd <- 0.2; lambda <- 1; eff.sd <- 0.3; eta.ka ~ 0.6 })
    model({ ka <- exp(tka + eta.ka); cl <- exp(tcl); v <- exp(tv); cp <- linCmt(); eff <- cp * exp(te0)
            cp ~ add(add.sd) + boxCox(lambda); eff ~ add(eff.sd) })
  }
  ui <- rxode2::rxode2(rmod)
  cl <- .rpemClassify(ui)
  expect_equal(cl$endpt$errType, c(3L, 0L))   # cp TBS (Box-Cox), eff additive
  rf <- .rpemFit(ui, dat, rpemControl(nGauss = 300L, nMH = 80000L, mhBurn = 8000L,
                                      niter = 30L, collect = 12L, seed = 300L, cores = 4L))
  # residuals recover the truth: transformed-scale add.sd, the Box-Cox lambda, and
  # the other endpoint's additive sd.
  expect_equal(rf$endptSd[1], addT, tolerance = 0.04)    # cp add.sd (transformed scale)
  expect_equal(rf$endptProp[1], lamT, tolerance = 0.2)   # cp lambda
  expect_equal(rf$endptSd[2], effT, tolerance = 0.06)    # eff add.sd
})
