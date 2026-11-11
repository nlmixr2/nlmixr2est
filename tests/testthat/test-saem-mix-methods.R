test_that("every saem mixProbMethod runs through the SA covariance phase (#1297)", {
  skip_on_cran()

  .testSeed(1001)
  nSub <- 40L
  grp <- sample.int(2L, nSub, TRUE, prob = c(0.45, 0.55))
  sim <- rxode2::rxode2({
    cl <- CLI
    d/dt(depot) <- -1.1 * depot
    d/dt(center) <- 1.1 * depot - cl / 20 * center
    cp <- center / 20
  })
  ev <- rxode2::et(rxode2::et(amt = 320, cmt = "depot"), c(0.5, 1, 2, 4, 8, 12, 24))
  obs <- do.call(
    rbind,
    lapply(seq_len(nSub), function(i) {
      s <- rxode2::rxSolve(sim, params = c(CLI = c(1, 8)[grp[i]]), ev, returnType = "data.frame")
      data.frame(ID = i, TIME = s$time, DV = s$cp + rnorm(nrow(s), 0, 0.05), AMT = 0, EVID = 0)
    })
  )
  dat <- rbind(data.frame(ID = seq_len(nSub), TIME = 0, DV = NA_real_, AMT = 320, EVID = 1), obs)
  dat <- dat[order(dat$ID, dat$TIME, -dat$EVID), ]

  mixmod <- function() {
    ini({
      tka <- fix(log(1.1)); tcl1 <- fix(log(1.0)); tcl2 <- fix(log(8.0))
      tv <- fix(log(20)); p1 <- 0.45; eta.cl ~ 0.01; add.sd <- fix(0.05)
    })
    model({
      ka <- exp(tka)
      cl <- mix(exp(tcl1 + eta.cl), p1, exp(tcl2 + eta.cl))
      v <- exp(tv)
      linCmt() ~ add(add.sd)
    })
  }

  for (.m in c("regress", "regularized", "annealed")) {
    fit <- suppressWarnings(nlmixr2(
      mixmod,
      dat,
      est = "saem",
      saemControl(
        print = 0,
        nBurn = 50L,
        nEm = 50L,
        covMethod = "sa",
        nSaCov = 20L,
        calcTables = FALSE,
        mixProbMethod = .m
      )
    ))
    .p <- unname(fit$env$mixProbabilities)
    expect_true(all(is.finite(.p)), info = .m)
    expect_equal(sum(.p), 1, tolerance = 1e-8, info = .m)
  }
})
