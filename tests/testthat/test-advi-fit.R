## est="advi" end-to-end (mean-field, point-estimate): the C++ optimization loop
## runs, the ELBO trend increases, and the population estimates are finite and
## physically sane.  Also checks same-seed reproducibility of the loop.

nmTest({
  test_that("est='advi' runs end-to-end and the ELBO increases (theo_sd)", {
    one.cmt <- function() {
      ini({ tka <- 0.45; tcl <- 1; tv <- 3.45; eta.ka ~ 0.6; add.sd <- 0.7 })
      model({ ka <- exp(tka + eta.ka); cl <- exp(tcl); v <- exp(tv)
        d/dt(depot) <- -ka*depot; d/dt(center) <- ka*depot - cl/v*center
        cp <- center/v; cp ~ add(add.sd) })
    }
    fit <- suppressMessages(suppressWarnings(
      nlmixr2(one.cmt, nlmixr2data::theo_sd, est = "advi",
              control = adviControl(iters = 150L, print = 0L))))
    expect_s3_class(fit, "nlmixr2advi")

    ## ELBO trend: mean of last decile beats mean of first decile
    e <- fit$elbo
    n <- length(e); d <- max(1L, n %/% 10L)
    expect_gt(mean(e[(n - d + 1L):n]), mean(e[1:d]))

    ## population estimates finite + physically sane (tka in a reasonable range)
    expect_true(all(is.finite(fit$theta)))
    expect_true(all(is.finite(fit$popOmega)) && all(fit$popOmega > 0))
    expect_gt(fit$theta[1], -1); expect_lt(fit$theta[1], 2)   # tka ~ log(ka)

    ## variational posterior params finite
    expect_equal(dim(fit$mu), c(length(unique(nlmixr2data::theo_sd$ID)), 1L))
    expect_true(all(is.finite(fit$mu)) && all(is.finite(fit$omega)))
  })

  test_that("est='advi' is reproducible under a fixed seed", {
    one.cmt <- function() {
      ini({ tka <- 0.45; tcl <- 1; tv <- 3.45; eta.ka ~ 0.6; add.sd <- 0.7 })
      model({ ka <- exp(tka + eta.ka); cl <- exp(tcl); v <- exp(tv)
        d/dt(depot) <- -ka*depot; d/dt(center) <- ka*depot - cl/v*center
        cp <- center/v; cp ~ add(add.sd) })
    }
    ctl <- adviControl(iters = 40L, seed = 7L, print = 0L)
    f1 <- suppressMessages(suppressWarnings(
      nlmixr2(one.cmt, nlmixr2data::theo_sd, est = "advi", control = ctl)))
    f2 <- suppressMessages(suppressWarnings(
      nlmixr2(one.cmt, nlmixr2data::theo_sd, est = "advi", control = ctl)))
    expect_identical(f1$theta, f2$theta)
    expect_identical(f1$elbo, f2$elbo)
    expect_identical(f1$mu, f2$mu)
  })
})
