nmTest({
  .emaxLbfgsb3c <- function() {
    ini({
      E0 <- 0.5
      Em <- 0.5
      E50 <- 2
      g <- fix(2)
    })
    model({
      v <- E0 + Em * time^g / (E50^g + time^g)
      ll(bin) ~ DV * v - log(1 + exp(v))
    })
  }

  .dsnLbfgsb3c <- function() {
    set.seed(42)
    dsn <- data.frame(i = 1:1000)
    dsn$time <- exp(rnorm(1000))
    dsn$DV <- rbinom(1000, 1, exp(-1 + dsn$time) / (1 + exp(-1 + dsn$time)))
    dsn
  }

  test_that("lbfgsb3cControl() keeps maxit", {
    expect_equal(lbfgsb3cControl()$maxit, 10000L)
    expect_equal(lbfgsb3cControl(maxit = 20)$maxit, 20L)
    expect_error(lbfgsb3cControl(maxit = 0))
    expect_equal(do.call(lbfgsb3cControl, lbfgsb3cControl(maxit = 7L))$maxit, 7L)
  })

  test_that("est='lbfgsb3c' runs in C++ without R objective callbacks", {
    skip_on_cran()
    .dsn <- .dsnLbfgsb3c()
    .ref <- suppressMessages(
      nlmixr2(.emaxLbfgsb3c, .dsn, est = "n1qn1", n1qn1Control(print = 0))
    )
    # The optimizer must never call back into R for fn/gr.
    local_mocked_bindings(
      .nlmixrOptimFunC = function(...) stop("R objective callback used"),
      .nlmixrOptimGradC = function(...) stop("R gradient callback used")
    )
    .ret <- suppressMessages(
      nlmixr2(.emaxLbfgsb3c, .dsn, est = "lbfgsb3c", lbfgsb3cControl(print = 0, returnLbfgsb3c = TRUE))
    )
    expect_equal(.ret$convergence, 0L)
    expect_match(.ret$message, "CONVERGENCE")
    expect_gt(.ret$counts[1], 0L)
    expect_true(is.finite(.ret$value))
    expect_named(.ret$grad, c("E0", "Em", "E50"))
    expect_named(.ret$par, c("E0", "Em", "E50"))

    # sigdig=3's factr stops ~1.6 OFV short here; sigdig=6 reaches the n1qn1 optimum
    .fit <- suppressMessages(
      nlmixr2(.emaxLbfgsb3c, .dsn, est = "lbfgsb3c", lbfgsb3cControl(print = 0, sigdig = 6))
    )
    expect_s3_class(.fit, "nlmixr2.lbfgsb3c")
    expect_equal(.fit$objf, .ref$objf, tolerance = 1e-4)
  })

  test_that("nlmLbfgsb3cFit() refuses an unloaded problem", {
    expect_error(nlmLbfgsb3cFit(1, -Inf, Inf, list()), "not loaded")
  })

  test_that("est='lbfgsb3c' honors maxit", {
    skip_on_cran()
    .ret <- suppressWarnings(suppressMessages(
      nlmixr2(
        .emaxLbfgsb3c,
        .dsnLbfgsb3c(),
        est = "lbfgsb3c",
        lbfgsb3cControl(print = 0, maxit = 2L, returnLbfgsb3c = TRUE)
      )
    ))
    expect_equal(.ret$convergence, 1L)
    expect_match(.ret$message, "Maximum number of iterations")
  })

  test_that("FOCEi outerOpt='lbfgsb3c' reports the L-BFGS-B exit message", {
    skip_on_cran()
    .one <- function() {
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
    .fit <- suppressWarnings(suppressMessages(
      nlmixr2(
        .one,
        nlmixr2data::theo_sd,
        est = "focei",
        foceiControl(outerOpt = "lbfgsb3c", fast = TRUE, print = 0, covMethod = "", calcTables = FALSE)
      )
    ))
    expect_true(is.finite(.fit$objf))
    expect_match(.fit$env$message, "CONVERGENCE|Maximum|WARNING|ABNORMAL")
  })
})
