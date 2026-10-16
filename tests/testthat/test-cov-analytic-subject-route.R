# The per-subject augmented solve, the fallback when the batched one fails
# (forced with FOCEI_NO_FD3_BATCH = 1), gives the batched route's covariance on
# each assembly: add/prop FOCEI, add + prop FOCE, and the general (f,R) assembly
# under FOCEI and FOCE, whose per-subject route has to rebuild the sigma slots
# from rsig as the batched one does.  Both are evaluated at one fit's estimates
# (maxOuterIterations = 0).  The routes difference the third-order tensor
# differently (Richardson against Shi); the largest measured difference in an SE
# is 1.1e-6, so 1e-5 is about ten times it.  On the (f,R) assembly the
# per-subject route used to error ("number of items to replace is not a multiple
# of replacement length").

nmTest({
  .srBase <- function() {
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
      d/dt(depot) <- -ka * depot
      d/dt(center) <- ka * depot - cl / v * center
      cp <- center / v
      cp ~ add(add.sd)
    })
  }

  .srAddProp <- function() {
    ini({
      tka <- 0.45
      tcl <- 1
      tv <- 3.45
      eta.ka ~ 0.6
      eta.cl ~ 0.3
      eta.v ~ 0.1
      add.sd <- 0.5
      prop.sd <- 0.1
    })
    model({
      ka <- exp(tka + eta.ka)
      cl <- exp(tcl + eta.cl)
      v <- exp(tv + eta.v)
      d/dt(depot) <- -ka * depot
      d/dt(center) <- ka * depot - cl / v * center
      cp <- center / v
      cp ~ add(add.sd) + prop(prop.sd)
    })
  }

  .srC1 <- function() {
    ini({
      tka <- 0.45
      tcl <- 1
      tv <- 3.45
      eta.ka ~ 0.6
      eta.cl ~ 0.3
      eta.v ~ 0.1
      add.sd <- 0.5
      prop.sd <- 0.1
    })
    model({
      ka <- exp(tka + eta.ka)
      cl <- exp(tcl + eta.cl)
      v <- exp(tv + eta.v)
      d/dt(depot) <- -ka * depot
      d/dt(center) <- ka * depot - cl / v * center
      cp <- center / v
      cp ~ add(add.sd) + prop(prop.sd) + combined1()
    })
  }

  .srFit <- function(model, est, ...) {
    suppressMessages(suppressWarnings(nlmixr2(
      model,
      nlmixr2data::theo_sd,
      est,
      foceiControl(print = 0L, calcTables = FALSE, ...),
      table = tableControl(cwres = FALSE, npde = FALSE)
    )))
  }

  .srCase <- list(
    addFocei = list(model = .srBase, est = "focei"),
    addPropFoce = list(model = .srAddProp, est = "foce"),
    generalFocei = list(model = .srC1, est = "focei"),
    generalFoce = list(model = .srC1, est = "foce")
  )

  test_that("the per-subject solve gives the batched solve's analytic covariance", {
    skip_on_cran()
    for (.label in names(.srCase)) {
      .c <- .srCase[[.label]]
      .fit <- .srFit(.c$model, .c$est, covMethod = "")
      .batched <- .srFit(.fit$finalUi, .c$est, covMethod = "analytic", maxOuterIterations = 0L)
      .subject <- local({
        withr::local_envvar(FOCEI_NO_FD3_BATCH = "1")
        .srFit(.fit$finalUi, .c$est, covMethod = "analytic", maxOuterIterations = 0L)
      })
      expect_false(exists(".analyticReason", envir = .batched$env, inherits = FALSE), label = .label)
      expect_false(exists(".analyticReason", envir = .subject$env, inherits = FALSE), label = .label)
      expect_identical(.subject$covMethod, .batched$covMethod, label = .label)
      expect_equal(sqrt(diag(.subject$cov)), sqrt(diag(.batched$cov)), tolerance = 1e-5, label = .label)
    }
  })
})
