nmTest({
  # #1286: a warm-started ETA can fall into a secondary mode during the outer
  # search and stay there, so the final point scored above an earlier
  # evaluation of the same theta.  Each subject's inner problem here has two
  # modes, eta = cc +/- sqrt(DV); the prior favors the one nearer 0.
  bi <- function() {
    ini({
      la <- 0
      cc <- 0.2
      addSd <- 0.1
      eta.x ~ 1
    })
    model({
      f <- exp(la) * (eta.x - cc)^2
      f ~ add(addSd)
    })
  }
  set.seed(1)
  n <- 20
  d <- data.frame(ID = rep(1:n, each = 4), TIME = rep(1:4, n))
  d$DV <- rep(runif(n, 0.8, 1.2), each = 4) + rnorm(nrow(d), 0, 0.1)

  # Evaluate the start, then a point that moves the ETAs past the second mode,
  # then the start again (warm, so they stay there); return `ret(par, par2)`.
  .jumpOpt <- function(ret) {
    function(par, fn, gr, lower = -Inf, upper = Inf, control = list(), ...) {
      par2 <- par
      par2[2] <- par2[2] + 0.5
      fn(par)
      fn(par2)
      fn(par)
      list(x = ret(par, par2), convergence = 0L, message = "jump")
    }
  }
  .fit <- function(opt, tol) {
    suppressMessages(nlmixr2(
      bi,
      d,
      "focei",
      control = foceiControl(
        print = 0,
        covMethod = "",
        calcTables = FALSE,
        outerOpt = opt,
        finalEtaTol = tol,
        resetEtaP = 1e-10
      )
    ))
  }

  test_that("the final ETAs are re-solved from the best evaluation's", {
    skip_on_cran()
    opt <- .jumpOpt(function(par, par2) par)
    f0 <- .fit(opt, Inf)
    f1 <- .fit(opt, 1)
    first <- f1$parHistData
    first <- first$objf[first$type == "Unscaled"][1]
    # without the re-solve the fit ends in the secondary mode
    expect_gt(f0$objf - first, 5)
    expect_equal(f0$env$finalEtaMoved, 0L)
    expect_gt(f1$env$finalEtaMoved, 0L)
    expect_equal(f1$objf, first, tolerance = 1e-3)
    expect_true(all(f1$eta$eta.x < 0))
    expect_true(any(grepl("re-solved from the best iterate", f1$runInfo)))
    expect_false(any(grepl("not converged", f1$runInfo)))
  })

  test_that("a final point above the best evaluated one is replaced by it", {
    skip_on_cran()
    f <- .fit(.jumpOpt(function(par, par2) par2), 1)
    expect_true(f$env$finalBestX)
    expect_equal(unname(fixef(f)["cc"]), 0.2, tolerance = 1e-6)
    expect_true(any(grepl("replaced by the lowest evaluated point", f$runInfo)))
  })

  test_that("finalEtaTol is validated", {
    expect_error(foceiControl(finalEtaTol = -1))
    expect_equal(foceiControl()$finalEtaTol, 1)
    expect_equal(foceiControl(finalEtaTol = Inf)$finalEtaTol, Inf)
  })
})
