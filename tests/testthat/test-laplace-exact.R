nmTest({
  test_that("laplaceControl() asks for the exact (log-density) Laplace Hessian", {
    expect_true(laplaceControl()$laplaceLlik)
    expect_true(mlaplaceControl()$laplaceLlik)
    expect_true(ilaplaceControl()$laplaceLlik)
    # the analytic conditional Hessian (flaplace) is the same objective by
    # another route, so it does not also force the log-density endpoints
    .c <- getValidNlmixrCtl(structure(list(NULL), class = "flaplace"))
    expect_false(.c$laplaceLlik)
    # a laplaceControl() handed to est="focei" keeps FOCEi: the flag is
    # laplace's own default, not a user setting, so it does not carry over
    .f <- getValidNlmixrCtl(structure(list(laplaceControl()), class = "focei"))
    expect_false(isTRUE(.f$laplaceLlik))
  })

  test_that("est='laplace' is NONMEM's LAPLACE INTERACTION objective (Wang 2007 Eq. 15)", {
    skip_on_cran()
    set.seed(42)
    tt <- c(0.5, 1, 2, 4, 8)
    d <- do.call(
      rbind,
      lapply(1:12, function(i) {
        a <- 10 * exp(rnorm(1, 0, sqrt(0.1)))
        k <- 0.3 * exp(rnorm(1, 0, sqrt(0.1)))
        f <- a * exp(-k * tt)
        data.frame(ID = i, TIME = tt, DV = f * (1 + 0.15 * rnorm(5)))
      })
    )
    th <- c(9, 0.25)
    om <- c(0.15, 0.2)
    pb <- 0.2
    # Wang (2007) Eq. 15, without the n*log(2*pi) constant: exact mode, exact
    # (numerical) second derivative of the full conditional -2LL in eta,
    # residual-variance terms included
    wang15 <- 0
    for (i in unique(d$ID)) {
      di <- d[d$ID == i, ]
      l <- function(eta) {
        f <- th[1] * exp(eta[1]) * exp(-th[2] * exp(eta[2]) * di$TIME)
        v <- (f * pb)^2
        sum(log(v) + (di$DV - f)^2 / v)
      }
      oi <- diag(1 / om)
      o <- optim(
        c(0, 0),
        function(e) l(e) + sum(e * (oi %*% e)),
        method = "BFGS",
        control = list(reltol = 1e-14, maxit = 1000)
      )
      e <- o$par
      h <- numDeriv::hessian(l, e)
      wang15 <- wang15 + l(e) + sum(log(om)) + sum(e * (oi %*% e)) + as.numeric(determinant(h / 2 + oi)$modulus)
    }

    model <- function() {
      ini({
        ta <- log(9)
        tk <- log(0.25)
        pb <- 0.2
        eta.a ~ 0.15
        eta.k ~ 0.2
      })
      model({
        a <- exp(ta + eta.a)
        k <- exp(tk + eta.k)
        f <- a * exp(-k * time)
        f ~ prop(pb)
      })
    }
    .ctl <- list(print = 0L, calcTables = FALSE, covMethod = "", maxOuterIterations = 0L)
    fitL <- .nlmixr(model, d, "laplace", do.call(laplaceControl, .ctl))
    fitF <- .nlmixr(model, d, "focei", do.call(foceiControl, .ctl))
    expect_equal(fitL$objf, wang15, tolerance = 1e-4)
    # FOCEi (Wang Eq. 16) is a different objective under proportional error
    expect_gt(abs(fitF$objf - wang15), 0.01)
  })
})
