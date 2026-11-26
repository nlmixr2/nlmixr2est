nmTest({
  test_that("an outer optimizer may ask for the gradient before the objective", {
    skip_on_cran()
    skip_if_not_installed("nlmixr2data")
    # the gradient's history row read the last objective iteration, which did not
    # exist yet (a segfault)
    grFirst <- function(par, fn, gr, ...) {
      g <- gr(par)
      list(par = par, value = fn(par), convergence = 0L, message = "gradient first", x = par, iterations = 0L)
    }
    m <- function() {
      ini({ tka <- 0.45; tcl <- 1; tv <- 3.45; eta.cl ~ 0.3; eta.v ~ 0.1; add.sd <- 0.7 })
      model({
        ka <- exp(tka); cl <- exp(tcl + eta.cl); v <- exp(tv + eta.v)
        d/dt(depot) <- -ka * depot; d/dt(center) <- ka * depot - cl / v * center
        cp <- center / v; cp ~ add(add.sd)
      })
    }
    for (fast in c(TRUE, FALSE)) {
      f <- suppressMessages(suppressWarnings(nlmixr2(
        m,
        nlmixr2data::theo_sd,
        "focei",
        foceiControl(print = 0L, covMethod = "", fast = fast, outerOpt = grFirst)
      )))
      expect_true(is.finite(f$objf))
      # the objective is evaluated first, so the first gradient is the initial one
      ph <- f$parHistData
      g <- ph[!(ph$type %in% c("Scaled", "Unscaled", "Back-Transformed")), ]
      expect_equal(g$iter[1], 1)
      if (!fast) expect_identical(as.character(g$type[1]), "Gill83 Gradient")
    }
  })
})
