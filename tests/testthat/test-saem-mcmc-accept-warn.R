nmTest({
  test_that("saemControl(mcmcAcceptWarn=) is validated", {
    expect_equal(saemControl()$mcmcAcceptWarn, c(0.05, 0.8))
    expect_equal(saemControl(mcmcAcceptWarn = NA)$mcmcAcceptWarn, NA)
    expect_equal(saemControl(mcmcAcceptWarn = NULL)$mcmcAcceptWarn, NA)
    for (.v in list(NA, c(0.2, 0.6))) {
      .d <- rxUiDeparse(saemControl(mcmcAcceptWarn = .v), "ctl")
      expect_equal(eval(.d[[3]])$mcmcAcceptWarn, .v)
    }
    .ctl <- saemControl(mcmcAcceptWarn = c(0.2, 0.6))
    expect_equal(do.call(saemControl, .ctl)$mcmcAcceptWarn, c(0.2, 0.6))
    expect_error(saemControl(mcmcAcceptWarn = c(0.8, 0.1)))
    expect_error(saemControl(mcmcAcceptWarn = c(-0.1, 0.8)))
    expect_error(saemControl(mcmcAcceptWarn = 0.1))
  })

  test_that("the note judges the tail of the estimation iterations", {
    .env <- new.env(parent = emptyenv())
    .env$saemControl <- saemControl(nBurn = 150, nEm = 100)
    .acc <- matrix(0.3, 300, 3, dimnames = list(NULL, c("prior", "rw", "coord")))
    .acc[, "prior"] <- 0.01 # kernel 1 is never judged
    .acc[1:150, "rw"] <- 0.01 # early-bad, late-good passes
    .acc[251:300, "coord"] <- 0.99 # rows past nBurn + nEm are the cov phase
    .env$mcmcAccept <- .acc
    expect_silent(.r <- .saemMcmcAcceptWarn(.env))
    expect_equal(.r, c(rw = 0.3, coord = 0.3))

    .env$mcmcAccept[151:250, "rw"] <- 0.02
    .env$mcmcAccept[151:250, "coord"] <- 0.95
    expect_warning(
      expect_warning(.saemMcmcAcceptWarn(.env), "MCMC acceptance low \\(rw 2%\\); try a smaller rmcmc"),
      "MCMC acceptance high \\(coord 95%\\); try a larger rmcmc"
    )
    .env$saemControl <- saemControl(nBurn = 150, nEm = 100, mcmcAcceptWarn = NA)
    expect_silent(.saemMcmcAcceptWarn(.env))
  })

  test_that("a fit notes a poorly tuned chain in $runInfo", {
    skip_on_cran()
    .mod <- function() {
      ini({
        tka <- 0.45
        tcl <- 1
        tv <- 3.45
        eta.ka ~ 0.6
        eta.cl ~ 0.3
        eta.v ~ 0.1
        add.sd <- 0.7
      })
      model({
        ka <- exp(tka + eta.ka)
        cl <- exp(tcl + eta.cl)
        v <- exp(tv + eta.v)
        linCmt() ~ add(add.sd)
      })
    }
    .fit <- function(r) {
      suppressMessages(nlmixr2(
        .mod(),
        theo_sd,
        est = "saem",
        control = saemControl(print = 0, nBurn = 50, nEm = 50, seed = 42L, covMethod = "", rmcmc = r)
      ))
    }
    .isAcc <- function(f) grepl("^MCMC acceptance", f$runInfo)
    .bad <- .fit(10)
    expect_true(any(grepl("^MCMC acceptance low .*try a smaller rmcmc$", .bad$runInfo)))
    expect_true(all(nchar(.bad$runInfo[.isAcc(.bad)]) < 75L))
    expect_false(any(.isAcc(.fit(0.5))))
  })
})
