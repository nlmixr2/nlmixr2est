# The post-fit covariance hooks of nlmixr2Est0() (R/nlmixr2Est.R).

nmTest({
  .oneCmt <- function() {
    ini({
      tka <- 0.45
      tcl <- 1.0
      tv <- 3.45
      add.sd <- 0.7
      eta.ka ~ 0.6
      eta.cl ~ 0.3
      eta.v ~ 0.1
    })
    model({
      ka <- exp(tka + eta.ka)
      cl <- exp(tcl + eta.cl)
      v <- exp(tv + eta.v)
      linCmt() ~ add(add.sd)
    })
  }

  .fitDeferred <- function(covMethod) {
    suppressMessages(nlmixr2(
      .oneCmt,
      nlmixr2data::theo_sd,
      est = "focei",
      control = foceiControl(print = 0, covMethod = covMethod, calcTables = FALSE)
    ))
  }

  test_that("a deferred sa/imp covariance that cannot be computed is reported", {
    local_mocked_bindings(.covRecompute = function(fit, method, control = NULL) NULL)
    expect_warning(
      .fit <- .fitDeferred("imp"),
      "\"imp\" covariance could not be computed; none installed",
      fixed = TRUE
    )
    expect_null(.fit$cov)
    expect_null(.fit$env$covOptions$imp)
  })

  test_that("a deferred sa/imp covariance that fell back to another one says so", {
    local_mocked_bindings(.covRecompute = function(fit, method, control = NULL) {
      .n <- c("tka", "tcl", "tv", "add.sd")
      list(cov = matrix(diag(0.01, 4), 4, 4, dimnames = list(.n, .n)), covMethod = "linFim")
    })
    expect_warning(
      .fit <- .fitDeferred("sa"),
      "\"sa\" covariance could not be computed; installed \"linFim\"",
      fixed = TRUE
    )
    expect_identical(.fit$covMethod, "linFim")
    expect_equal(unname(sqrt(diag(.fit$cov))), rep(0.1, 4))
    # the options are those of the covariance requested, not the one installed
    expect_null(.fit$env$covOptions$sa)
  })
})
