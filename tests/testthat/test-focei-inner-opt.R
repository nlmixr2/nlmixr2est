# innerOpt="steihaug" (truncated-CG trust region) and innerOpt="newuoa"
# (derivative-free) solve the same per-subject eta problem as "trust".

nmTest({
  # nolint start: object_usage_linter.
  .mod <- function() {
    ini({
      tka <- 0.45
      tcl <- 1
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
  # nolint end
  .etaCols <- c("eta.ka", "eta.cl", "eta.v")
  .fit <- function(innerOpt, maxOuter = 0L, ...) {
    .nlmixr(.mod, nlmixr2data::theo_sd, est = "focei",
            control = foceiControl(print = 0, innerOpt = innerOpt,
                                   maxOuterIterations = maxOuter,
                                   covMethod = "", calcTables = FALSE, ...))
  }

  test_that("steihaug and newuoa reach trust's posthoc etas and objective", {
    .tr <- .fit("trust", trustPolish = TRUE)
    for (.o in c("steihaug", "newuoa")) {
      for (.p in c(FALSE, TRUE)) {
        .args <- list(.o)
        .args[[paste0(.o, "Polish")]] <- .p
        .f <- do.call(.fit, .args)
        expect_true(.f$env$nTrustInner[["calls"]] > 0, info = .o)
        expect_equal(.f$objf, .tr$objf, tolerance = 1e-3, info = paste(.o, .p))
        expect_equal(as.matrix(.f$eta[, .etaCols]), as.matrix(.tr$eta[, .etaCols]),
                     tolerance = 1e-2, ignore_attr = TRUE, info = paste(.o, .p))
      }
    }
  })

  test_that("steihaug fits like trust", {
    .tr <- .fit("trust", maxOuter = 100L)
    .st <- .fit("steihaug", maxOuter = 100L, steihaugPolish = TRUE)
    expect_equal(.st$objf, .tr$objf, tolerance = 0.05)
  })

  test_that("inner optimizer options validate and round-trip", {
    expect_error(foceiControl(innerOpt = "newuoa", hessianMethod = "bfgs"),
                 "hessianMethod")
    expect_error(foceiControl(innerOpt = "steihaug", hessianMethod = "bfgs"), NA)
    .ctl <- foceiControl(innerOpt = "steihaug", steihaugPolish = TRUE,
                         newuoaRhoend = 1e-4)
    expect_equal(rxUiDeparse.foceiControl(.ctl, "var"),
                 quote(var <- foceiControl(innerOpt = "steihaug",
                                           steihaugPolish = TRUE,
                                           newuoaRhoend = 1e-4)))
    # the default rhoend follows the ODE tolerance
    expect_equal(foceiControl()$newuoaRhoend, sqrt(max(foceiControl()$rxControl$rtol)))
  })
})
