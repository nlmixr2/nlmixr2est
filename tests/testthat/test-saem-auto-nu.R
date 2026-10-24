test_that("a residual error model is traced back to its etas", {
  .ui <- function(extraIni, lines) {
    eval(parse(
      text = sprintf(
        "function() {
        ini({ tka <- 0.45; tcl <- 1; tv <- 3.45; add.sd <- 0.5; eta.ka ~ 0.2; %s })
        model({ ka <- exp(tka + eta.ka); cl <- exp(tcl); v <- exp(tv); cp <- linCmt(); %s })
      }",
        extraIni,
        lines
      )
    ))()
  }
  .has <- function(ui) .saemModeledResidHasEta(ui, .saemModeledResidualCond(ui))

  expect_true(.has(.ui("eta.sd ~ 0.1", "a <- add.sd * exp(eta.sd); cp ~ add(a)")))
  # through intermediate variables
  expect_true(.has(.ui("eta.sd ~ 0.1", "l <- log(add.sd) + eta.sd; a <- exp(l); cp ~ add(a)")))
  # an explicit + dnorm() is the same likelihood
  .ui2 <- .ui("eta.sd ~ 0.1", "a <- add.sd * exp(eta.sd); cp ~ add(a) + dnorm()")
  expect_true(.saemModeledResidHasEta(.ui2, .saemModeledResidualCond(.ui2, c("norm", "dnorm"))))
  # a covariate alone does not slow the chains
  expect_false(.has(.ui("wt.sd <- 0.01", "a <- add.sd + wt.sd * WT; cp ~ add(a)")))
})

test_that("only an unset saem nu is raised", {
  expect_equal(.saemAutoNu(NULL)$mcmc$nu, c(4, 4, 4))
  expect_equal(.saemAutoNu(saemControl(print = 0))$mcmc$nu, c(4, 4, 4))
  expect_null(.saemAutoNu(saemControl(nu = c(6, 1, 1))))
  expect_null(.saemAutoNu(saemControl(nu = c(2, 2, 2))))
  expect_null(.saemAutoNu(list(print = 0)))
  # saem rebuilds its control through do.call(saemControl, ...)
  .rebuilt <- do.call(saemControl, unclass(saemControl()))
  expect_true(.rebuilt$mcmc$nuAuto)
  expect_false(do.call(saemControl, unclass(saemControl(nu = c(2, 2, 2))))$mcmc$nuAuto)
})
