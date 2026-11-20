nmTest({
  test_that("saemControl(rmcmc=) is validated and survives a rebuild", {
    expect_equal(saemControl()$mcmc$rmcmc, 0.5)
    .ctl <- saemControl(rmcmc = 2)
    expect_equal(.ctl$mcmc$rmcmc, 2)
    # saem rebuilds its control through do.call(saemControl, ...)
    expect_equal(do.call(saemControl, .ctl)$mcmc$rmcmc, 2)
    expect_error(saemControl(rmcmc = -1))
    expect_error(saemControl(rmcmc = c(1, 2)))
    expect_error(saemControl(rmcmc = Inf))
    # zero proposals never move the chain
    expect_error(saemControl(rmcmc = 0), "positive")
    .d <- rxUiDeparse(saemControl(rmcmc = 0.2, nmc = 5), "ctl")
    expect_equal(eval(.d[[3]])$mcmc$rmcmc, 0.2)
    expect_equal(eval(rxUiDeparse(saemControl(), "ctl")[[3]])$mcmc$rmcmc, 0.5)
  })

  test_that("rmcmc reaches the sampler", {
    skip_on_cran()
    .mod <- function() {
      ini({
        tka <- 0.45; tcl <- 1; tv <- 3.45
        eta.ka ~ 0.6; eta.cl ~ 0.3; eta.v ~ 0.1
        add.sd <- 0.7
      })
      model({
        ka <- exp(tka + eta.ka); cl <- exp(tcl + eta.cl); v <- exp(tv + eta.v)
        linCmt() ~ add(add.sd)
      })
    }
    .fit <- function(r) {
      suppressMessages(suppressWarnings(nlmixr2(
        .mod(),
        theo_sd,
        est = "saem",
        control = saemControl(print = 0, nBurn = 20, nEm = 20, seed = 42L, covMethod = "", rmcmc = r)
      )))
    }
    .lo <- .fit(0.1)
    .hi <- .fit(5)
    # wider random-walk proposals are rejected more often
    expect_gt(mean(.lo$mcmcAccept[, "rw"]), mean(.hi$mcmcAccept[, "rw"]))
    expect_gt(mean(.lo$mcmcAccept[, "coord"]), mean(.hi$mcmcAccept[, "coord"]))
  })
})
