## The FOCEi objective must be a function of `theta` alone.
##
## The inner eta reset / eta nudge machinery carries per-subject state (the
## warm-started EBEs) across outer iterations.  When that state is allowed to
## replace -- rather than compete with -- a converged EBE, the same `theta` can
## return values far apart, which corrupts the outer optimizer's model and makes
## it stall, oscillate, and exit at a point worse than one already visited.
##
## The gate below is the invariant itself, so it cannot rot into prose: the
## objective a fit reports must not be worse than a cold-start (etas from zero)
## evaluation at the fit's own final estimates.

nmTest({
  test_that("FOCEi objective does not depend on the inner eta history", {

    ## ---- simulate: 1-cmt oral, Michaelis-Menten elimination, 4 etas ---------
    ## Ingredients that make the reset fire often:
    ##  - four random effects,
    ##  - true BSV (0.6) above diagOmegaBoundUpper * initial omega (5 * 0.1),
    ##    so omega pins at its bound and standardized etas are inflated,
    ##  - sparse sampling with an LLOQ, so some subjects carry little eta
    ##    information and the inner optimizer stalls at its starting point.
    simMod <- rxode2::rxode2({
      ka   <- exp(lka   + eta.ka)
      vc   <- exp(lvc   + eta.vc)
      vmax <- exp(lvmax + eta.vmax)
      km   <- exp(lkm   + eta.km)
      cp   <- centr / vc
      d/dt(depot) <- -ka * depot
      d/dt(centr) <-  ka * depot - vmax * cp / (km + cp)
    })

    nsub  <- 32L
    doses <- rep(c(10, 50, 200), length.out = nsub)
    tobs  <- c(0.5, 1, 2, 4, 8, 12, 24, 48)
    ev <- do.call(rbind, lapply(seq_len(nsub), function(i) {
      rbind(data.frame(id = i, time = 0, amt = doses[i], evid = 1),
            data.frame(id = i, time = tobs, amt = 0, evid = 0))
    }))
    omTrue <- lotri::lotri(eta.ka + eta.vc + eta.vmax + eta.km ~
                             c(0.6, 0, 0.6, 0, 0, 0.6, 0, 0, 0, 0.6))
    sim <- rxode2::rxSolve(simMod, ev,
                           params = c(lka = log(0.8), lvc = log(30),
                                      lvmax = log(15), lkm = log(2)),
                           omega = omTrue, returnType = "data.frame",
                           addDosing = TRUE, seed = 20260727)

    ## rxSolve(addDosing=TRUE) marks observation rows evid = 2
    d <- sim[, c("id", "time", "amt", "evid", "cp")]
    names(d) <- c("ID", "TIME", "AMT", "EVID", "DV")
    d$EVID <- ifelse(d$EVID == 1, 1L, 0L)
    withr::with_seed(11, {
      n <- nrow(d)
      d$DV <- d$DV * (1 + rnorm(n, 0, 0.2)) + rnorm(n, 0, 0.05)
    })
    d$DV[d$EVID == 1] <- NA_real_
    LLOQ <- 0.05
    d$CENS <- ifelse(d$EVID == 0 & !is.na(d$DV) & d$DV < LLOQ, 1L, 0L)
    d$CENS[d$EVID == 1] <- NA_integer_
    d$DV[which(d$CENS == 1)] <- LLOQ

    expect_true(sum(d$CENS == 1, na.rm = TRUE) > 0)

    fitMod <- function() {
      ini({
        lka   <- log(0.4)
        lvc   <- log(50)
        lvmax <- log(8)
        lkm   <- log(1)
        eta.ka   ~ 0.1
        eta.vc   ~ 0.1
        eta.vmax ~ 0.1
        eta.km   ~ 0.1
        propSd <- 0.2
        addSd  <- 0.05
      })
      model({
        ka   <- exp(lka   + eta.ka)
        vc   <- exp(lvc   + eta.vc)
        vmax <- exp(lvmax + eta.vmax)
        km   <- exp(lkm   + eta.km)
        cp   <- centr / vc
        d/dt(depot) <- -ka * depot
        d/dt(centr) <-  ka * depot - vmax * cp / (km + cp)
        cp ~ prop(propSd) + add(addSd)
      })
    }

    fit <- suppressWarnings(
      nlmixr2(fitMod, d, est = "focei", control = list(print = 0L)))

    ## (1) the optimizer must not return a point worse than one it evaluated
    trace <- fit$parHistData$objf[fit$parHistData$type == "Unscaled"]
    expect_lt(fit$objDf$OBJF[1] - min(trace, na.rm = TRUE), 1)

    ## (2) the reported objective must not be inflated relative to a cold-start
    ##     evaluation at the same theta -- i.e. it is a function of theta alone.
    pf <- fit$parFixedDf
    th <- setNames(pf$Est, rownames(pf))
    om <- diag(fit$omega)
    om <- om[om > 0]
    cold <- suppressWarnings(nlmixr2(
      fitMod |> rxode2::ini(lka = th[["lka"]], lvc = th[["lvc"]],
                            lvmax = th[["lvmax"]], lkm = th[["lkm"]],
                            propSd = th[["propSd"]], addSd = th[["addSd"]],
                            eta.ka = om[[1]], eta.vc = om[[2]],
                            eta.vmax = om[[3]], eta.km = om[[4]]),
      d, est = "focei",
      control = list(print = 0L, maxOuterIterations = 0L,
                     covMethod = "", calcTables = FALSE)))

    expect_lt(fit$objDf$OBJF[1] - cold$objDf$OBJF[1], 1)
  })
})
