nmTest({

  ## FOCEi's mixture-proportion gradient is analytic: mixGrad() (src/inner.cpp)
  ## short-circuits the finite difference in numericGrad(), so a wrong value is
  ## invisible under the default derivative-free outer optimizer (bobyqa) and
  ## silently wrong under every gradient-based one.
  ##
  ## The exact gradient is
  ##   d(-2*log-lik)/d(mlogit theta_l) = -2 * sum_i (r_il - pi_l),
  ## whose stationary point is the EM condition pi_l == mean_i r_il.  The old
  ## form summed d(lik)/d(pi_l) and chained it through the DIAGONAL of the
  ## mexpit Jacobian, dropping both the off-diagonal terms (wrong from nMix == 3
  ## up) and the -2 objective scale (wrong at every nMix, sign included).  The
  ## flipped sign sent lbfgsb3c's first step uphill; at nMix >= 3 the wrong
  ## Jacobian compounds that and the proportions never moved off their initial
  ## values at all.  At nMix == 2 the old value is exactly the right one divided
  ## by -2 (the diagonal-only chain reduces to the correct expression there).

  .simMixData <- function(clTrue, pTrue, nSub, seed = 42L) {
    set.seed(seed)
    .grp <- sample.int(length(clTrue), nSub, TRUE, prob = pTrue)
    .sim <- rxode2::rxode2({
      ka <- 1.1
      cl <- CLI
      v <- 20
      d / dt(depot) <- -ka * depot
      d / dt(center) <- ka * depot - cl / v * center
      cp <- center / v
    })
    .ev <- rxode2::et(amt = 320, cmt = "depot")
    .ev <- rxode2::et(.ev, c(0.25, 0.5, 1, 2, 4, 8, 12, 24))
    .obs <- do.call(rbind, lapply(seq_len(nSub), function(i) {
      .s <- rxode2::rxSolve(.sim, params = c(CLI = clTrue[.grp[i]]), .ev,
                            returnType = "data.frame")
      data.frame(ID = i, TIME = .s$time,
                 DV = .s$cp + stats::rnorm(nrow(.s), 0, 0.05),
                 AMT = 0, EVID = 0)
    }))
    .dose <- data.frame(ID = seq_len(nSub), TIME = 0, DV = NA_real_,
                        AMT = 320, EVID = 1)
    .d <- rbind(.dose, .obs)
    list(data = .d[order(.d$ID, .d$TIME, -.d$EVID), ],
         emp = as.numeric(table(.grp)) / nSub)
  }

  ## Everything but the mixture proportions and one component's typical value
  ## is fix()ed, so mixGrad() drives essentially all of the outer optimization
  ## and a wrong value cannot be masked by the rest of the model improving the
  ## objective on its own.  It is deliberately NOT the only free parameter:
  ## with npars == 1 the outer optimization never reaches numericGrad(), so
  ## mixGrad() would not run at all and the test would be vacuous.
  .ctlGrad <- function() {
    foceiControl(print = 0, outerOpt = "lbfgsb3c", maxOuterIterations = 300L,
                 maxInnerIterations = 100L, covMethod = "", calcTables = FALSE)
  }
  .ctlNone <- function() {
    foceiControl(print = 0, maxOuterIterations = 0L,
                 maxInnerIterations = 100L, covMethod = "", calcTables = FALSE)
  }

  .checkMixGrad <- function(mod, dat, emp) {
    .f0 <- suppressWarnings(nlmixr2(mod, dat$data, "focei", .ctlNone()))
    .f <- suppressWarnings(nlmixr2(mod, dat$data, "focei", .ctlGrad()))
    .pi <- .f$env$mixProbabilities
    .rbar <- vapply(.f$env$mixList, function(z) mean(z$prob), numeric(1))
    ## it moved at all -- the sign check.  Under the old gradient the objective
    ## was bit-identical to the 0-iteration one.
    expect_true(.f$objective < .f0$objective - 1)
    ## it landed on the EM fixed point -- the Jacobian check
    expect_equal(unname(.pi), unname(.rbar), tolerance = 1e-3)
    ## and that fixed point is the truth the data was simulated from
    expect_equal(unname(.pi), emp, tolerance = 0.02)
  }

  test_that("focei estimates 2-component mixture proportions from the analytic gradient", {
    .dat <- .simMixData(c(1.0, 6.0), c(0.6, 0.4), 40L)
    .mod <- function() {
      ini({
        tka <- fix(log(1.1))
        tcl1 <- fix(log(1.0))
        tcl2 <- log(6.0)
        tv <- fix(log(20))
        p1 <- 0.25
        eta.cl ~ fix(0.01)
        add.sd <- fix(0.05)
      })
      model({
        ka <- exp(tka)
        cl <- mix(exp(tcl1 + eta.cl), p1, exp(tcl2 + eta.cl))
        v <- exp(tv)
        linCmt() ~ add(add.sd)
      })
    }
    .checkMixGrad(.mod, .dat, .dat$emp)
  })

  test_that("focei estimates 3-component mixture proportions from the analytic gradient", {
    ## nMix >= 3 is where the dropped off-diagonal mexpit Jacobian terms bite;
    ## at nMix == 2 the diagonal-only chain happens to reduce to the right
    ## expression (up to the -2 the old code also dropped).
    .dat <- .simMixData(c(1.0, 4.0, 12.0), c(0.5, 0.3, 0.2), 60L)
    .mod <- function() {
      ini({
        tka <- fix(log(1.1))
        tcl1 <- fix(log(1.0))
        tcl2 <- fix(log(4.0))
        tcl3 <- log(12.0)
        tv <- fix(log(20))
        p1 <- 0.2
        p2 <- 0.2
        eta.cl ~ fix(0.01)
        add.sd <- fix(0.05)
      })
      model({
        ka <- exp(tka)
        cl <- mix(exp(tcl1 + eta.cl), p1, exp(tcl2 + eta.cl), p2,
                  exp(tcl3 + eta.cl))
        v <- exp(tv)
        linCmt() ~ add(add.sd)
      })
    }
    .checkMixGrad(.mod, .dat, .dat$emp)
  })

})
