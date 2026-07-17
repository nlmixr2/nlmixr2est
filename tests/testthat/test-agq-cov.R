# Analytic AGQ covariance (covType="analytic", nAGQ > 1): the node terms, and the scope
# gates.  The gates matter because an ungated model fails SILENTLY -- it returns the nAGQ=1
# Laplace covariance stamped covMethod="analytic", with no error.

nmTest({
  .agq_cov_mod <- function() {
    ini({
      tka <- log(1.5); tcl <- log(2.7); tv <- log(31.5)
      eta.ka ~ 0.6; eta.cl ~ 0.3; eta.v ~ 0.1
      add.sd <- 0.7
    })
    model({
      ka <- exp(tka + eta.ka); cl <- exp(tcl + eta.cl); v <- exp(tv + eta.v)
      d/dt(depot)  <- -ka * depot
      d/dt(center) <-  ka * depot - cl / v * center
      cp <- center / v
      cp ~ add(add.sd)
    })
  }
  .agq_cov_1eta <- function() {          # neta == 1: the commonest AGQ case
    ini({
      tka <- log(1.5); tcl <- log(2.7); tv <- log(31.5)
      eta.cl ~ 0.3
      add.sd <- 0.7
    })
    model({
      ka <- exp(tka); cl <- exp(tcl + eta.cl); v <- exp(tv)
      d/dt(depot)  <- -ka * depot
      d/dt(center) <-  ka * depot - cl / v * center
      cp <- center / v
      cp ~ add(add.sd)
    })
  }

  test_that("the AGQ branch with a single node reproduces the FOCEI branch exactly", {
    skip_on_cran()
    skip_on_ci()
    fit <- suppressMessages(nlmixr(.agq_cov_mod, nlmixr2data::theo_sd, "focei",
                                   foceiControl(print = 0L, covMethod = "")))
    ui <- rxode2::rxUiDecompress(fit$finalUi)
    dirs <- paste0("ETA_", 1:3, "_")
    am <- .foceiAnalyticAugModelDirs(ui, dirs)
    skip_if(is.null(am))
    ef <- .foceiAnalyticErrFull(ui)
    tol <- .foceiAnalyticSolveTol(ui)
    th <- fit$theta
    thv <- setNames(as.numeric(th), paste0("THETA_", seq_along(th), "_"))
    ehat <- as.numeric(fit$eta[1, -1])
    s1 <- fit$dataSav[fit$dataSav$ID == 1, , drop = FALSE]
    obs <- s1[s1$EVID == 0, , drop = FALSE]
    E <- .foceiAnalyticSolveSubjectFD3(am, c(thv, setNames(ehat, paste0("ETA_", 1:3, "_"))),
                                       s1, obs$TIME, tol = tol)
    skip_if(is.null(E))
    E$y <- .foceiAnalyticTbsY(obs$DV, E$trans)
    Om <- fit$omega
    args <- list(E = E, ehat = ehat, Om = Om, ef = ef, neta = 3L, nth = 3L, nsg = 1L,
                 sgVar = ef$sgVar, omd = .omegaVarCovDeriv(Om, .foceiOmegaPairs(Om, ui$iniDf)),
                 ndir = 3L, dirTh = 1:3, Oi = solve(Om))
    Rfocei <- do.call(.foceiAnalyticSubjectR, args)
    # Single node at x=0, so eta_1 == ehat; hand back `E` to test the formula, not the solver.
    # With etaP = -H^-1 M the AGQ data half folds to d2Phi - M_a'H^-1 M_b (= FOCEI), leaving
    # Ragq - Rfocei == g'eta2 with g = Phi_eta(etahat).  Both are exact for the objective as
    # defined (Phi + 0.5*log|Ht| AT THE MODE, where g == 0); the tolerance covers this fit's
    # stationarity residual (measured max|g| ~ 3.5e-5, and `epsilon` does not shrink it).
    # A structural break in the collapse is O(1) and is caught comfortably.
    g1 <- .agq(3L, 1L)
    Ragq <- do.call(.foceiAnalyticSubjectR,
                    c(args, list(qx = g1$x, qw = g1$w, solveNode = function(etak) E)))
    expect_false(is.null(Ragq))
    expect_equal(Ragq, Rfocei, tolerance = 1e-4)
  })

  test_that("covType='analytic' engages for nAGQ>1 and for a single eta", {
    skip_on_cran()
    skip_on_ci()
    for (.n in c(2L, 3L)) {
      fit <- suppressMessages(nlmixr(.agq_cov_mod, nlmixr2data::theo_sd, "agq",
                                     agqControl(print = 0L, covMethod = "", nAGQ = .n)))
      r <- foceiCovAnalytic(fit)
      expect_false(is.null(r))
      expect_identical(r$method, "analytic")
      expect_true(all(is.finite(r$se)))
    }
    # neta == 1 was dead code: vapply(, numeric(1)) drops dim, vK[[k]][, p] threw, and the
    # tryCatch silently fell back to FD.  This must return a real analytic covariance.
    fit1 <- suppressMessages(nlmixr(.agq_cov_1eta, nlmixr2data::theo_sd, "agq",
                                    agqControl(print = 0L, covMethod = "", nAGQ = 3L)))
    r1 <- foceiCovAnalytic(fit1)
    expect_false(is.null(r1))
    expect_true(all(is.finite(r1$se)))
  })

  test_that("AGQ declines (message + NULL) for every scope the node terms do not cover", {
    skip_on_cran()
    skip_on_ci()
    # interaction=0 routes to the FOCE assembler, which takes no qx/qw and would silently
    # return the nAGQ=1 FOCE covariance for an AGQ fit.
    fit <- suppressMessages(nlmixr(.agq_cov_mod, nlmixr2data::theo_sd, "agq",
                                   agqControl(print = 0L, covMethod = "", nAGQ = 2L,
                                              interaction = 0L)))
    expect_message(r <- foceiCovAnalytic(fit), "adaptive Gaussian quadrature")
    expect_null(r)
  })
})
