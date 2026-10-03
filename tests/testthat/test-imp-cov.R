nmTest({
  .impCovModel <- function() {
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

  test_that("llikObs is that of the estimates, not of the last importance-sampling covariance leg", {
    # impComputeCov() scores every subject at its fixed importance samples, at
    # perturbed parameters, after the final MAP pass set llikObs at the estimates
    for (.est in c("impmap", "imp")) {
      .ctl <- function(covMethod) {
        impmapControl(print = 0L, nIter = 5L, isample = 100L, covMethod = covMethod, calcTables = FALSE)
      }
      .imp <- .nlmixr(.impCovModel, theo_sd, .est, .ctl("imp"))
      .none <- .nlmixr(.impCovModel, theo_sd, .est, .ctl(""))
      expect_true(is.environment(.imp$env) && exists("impCovThetaN", envir = .imp$env))
      expect_identical(.imp$llikObs, .none$llikObs, label = .est)
    }
  })
})
