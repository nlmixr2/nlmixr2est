## Regression test for issue 358: a bioavailability/lag/rate/dur property on the
## 'central' compartment of a first-order absorption linCmt() model is silently
## ignored while solving (the dose enters the implied 'depot' compartment).
## nlmixr2est now warns the user to move the property to 'depot'.

nmTest({
  test_that("issue 358: warn about dosing properties on the wrong linCmt() compartment", {
    .hook <- get(".preProcessLinCmtDosing", envir = asNamespace("nlmixr2est"))
    .warns <- function(fun) {
      .ui <- rxode2::rxode2(fun)
      tryCatch({
        .hook(.ui, "focei", NULL, NULL)
        FALSE
      }, warning = function(cnd) TRUE)
    }

    firstOrderFCentral <- function() {
      ini({tka <- .5; tcl <- -3.2; tv <- -1; beta_f <- -0.5; eta.ka ~ 1; add.err <- 0.1})
      model({
        ka <- exp(tka + eta.ka); cl <- exp(tcl); v <- exp(tv)
        f(central) = 1 + bio_cov * beta_f
        linCmt() ~ add(add.err)
      })
    }
    firstOrderAlagCentral <- function() {
      ini({tka <- .5; tcl <- -3.2; tv <- -1; lagC <- 0.5; eta.ka ~ 1; add.err <- 0.1})
      model({
        ka <- exp(tka + eta.ka); cl <- exp(tcl); v <- exp(tv)
        alag(central) = lagC
        linCmt() ~ add(add.err)
      })
    }
    firstOrderFDepot <- function() {
      ini({tka <- .5; tcl <- -3.2; tv <- -1; beta_f <- -0.5; eta.ka ~ 1; add.err <- 0.1})
      model({
        ka <- exp(tka + eta.ka); cl <- exp(tcl); v <- exp(tv)
        f(depot) = 1 + bio_cov * beta_f
        linCmt() ~ add(add.err)
      })
    }
    ivFCentral <- function() {
      ini({tcl <- -3.2; tv <- -1; beta_f <- -0.5; eta.cl ~ 1; add.err <- 0.1})
      model({
        cl <- exp(tcl + eta.cl); v <- exp(tv)
        f(central) = 1 + bio_cov * beta_f
        linCmt() ~ add(add.err)
      })
    }

    # first-order absorption: property on 'central' is ignored -> warn
    expect_true(.warns(firstOrderFCentral))
    expect_true(.warns(firstOrderAlagCentral))
    # property on the correct dosing compartment -> no warning
    expect_false(.warns(firstOrderFDepot))
    # IV linCmt doses 'central' directly, so f(central) is legitimate -> no warning
    expect_false(.warns(ivFCentral))
  })

  test_that("issue 358: the warning reaches the fit runInfo", {
    skip_on_cran()

    data <- nlmixr2data::theo_sd
    data$bio_cov <- as.integer(data$ID > 6)

    fCentral <- function() {
      ini({tka <- .5; tcl <- 1; tv <- 3.45; beta_f <- -0.5; eta.ka ~ 1; add.err <- 0.7})
      model({
        ka <- exp(tka + eta.ka); cl <- exp(tcl); v <- exp(tv)
        f(central) = 1 + bio_cov * beta_f
        linCmt() ~ add(add.err)
      })
    }

    fit <- suppressWarnings(nlmixr2(
      fCentral, data, est = "focei",
      control = foceiControl(print = 0, covMethod = "", calcTables = FALSE,
                             maxOuterIterations = 0, maxInnerIterations = 0)))

    expect_true(any(grepl("linCmt.*depot|depot.*f\\(central\\)", fit$runInfo)))
  })
})
