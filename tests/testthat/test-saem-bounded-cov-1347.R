nmTest({
  # #1347: the bounded transform hid a bounded covariate coefficient from saem's
  # mu-referenced covariate regression, so it was mis-estimated (and reported on
  # the internal scale).  saem now keeps it in the closed-form M-step and holds
  # it inside its bounds there.
  .d <- nlmixr2data::theo_sd
  .d$lwt <- log(.d$WT / 70)

  .base <- function() {
    ini({
      tka <- 0.45
      tcl <- 1
      tv <- 3.45
      wt.cl <- 0.5
      eta.ka ~ 0.6
      eta.cl ~ 0.3
      eta.v ~ 0.1
      add.sd <- 0.7
    })
    model({
      ka <- exp(tka + eta.ka)
      cl <- exp(tcl + wt.cl * lwt + eta.cl)
      v <- exp(tv + eta.v)
      linCmt() ~ add(add.sd)
    })
  }
  # wt.cl with the bounds lo/hi (and an estimate between them)
  .mod <- function(lo = -Inf, est = 0.5, hi = Inf) {
    .ui <- rxode2::rxode2(.base)
    .iniDf <- .ui$iniDf
    .w <- .iniDf$name == "wt.cl"
    .iniDf$lower[.w] <- lo
    .iniDf$est[.w] <- est
    .iniDf$upper[.w] <- hi
    rxode2::ini(.ui) <- .iniDf
    .ui
  }

  .ctl <- saemControl(nBurn = 50, nEm = 50, seed = 42, print = 0, covMethod = "", calcTables = FALSE)

  .fit <- function(...) {
    suppressMessages(suppressWarnings(nlmixr2(.mod(...), .d, est = "saem", control = .ctl)))
  }

  test_that("phi1 coefficient bounds follow the kernel's Plambda1 order", {
    .ui <- .mod(-5, 0.5, 5) |>
      model(v <- exp(tv + wt.v * WT + eta.v)) |>
      ini(wt.v = c(0, 0.1, 2))
    expect_setequal(.saemLambda1BoundNames(.ui), c("wt.cl", "wt.v"))
    .cfg <- list(i1 = 0:2, nlambda1 = 5L)
    .b <- .saemLambda1Bounds(.ui, .cfg)
    # columns ka, cl, v: each its intercept, then its covariates (lwt, WT)
    expect_equal(.b$lower, c(-Inf, -Inf, -5, -Inf, 0))
    expect_equal(.b$upper, c(Inf, Inf, 5, Inf, 2))
  })

  test_that("a bounded covariate coefficient is estimated as an unbounded one (#1347)", {
    skip_on_cran()
    .free <- .fit()
    .bnd <- .fit(-5, 0.5, 5)
    expect_equal(.testBoundedTransform(), c(pre = FALSE, post = FALSE))
    expect_equal(.bnd$theta[["wt.cl"]], .free$theta[["wt.cl"]], tolerance = 1e-6)
    expect_equal(unname(.bnd$parFixedDf[["Estimate"]][rownames(.bnd$parFixedDf) == "wt.cl"]), .bnd$theta[["wt.cl"]])
    expect_equal(
      unname(.bnd$parFixedDf[["Back-transformed"]][rownames(.bnd$parFixedDf) == "wt.cl"]),
      .bnd$theta[["wt.cl"]]
    )
    expect_false(any(grepl("wt.cl", .bnd$runInfo, fixed = TRUE) & grepl("ignored", .bnd$runInfo)))
  })

  test_that("a binding bound holds the covariate coefficient at the bound (#1347)", {
    skip_on_cran()
    .free <- .fit()
    .hi <- .free$theta[["wt.cl"]] - 0.25
    .bnd <- .fit(-5, .hi - 0.1, .hi)
    expect_equal(.bnd$theta[["wt.cl"]], .hi)
    .lo <- .free$theta[["wt.cl"]] + 0.25
    .bnd <- .fit(.lo, .lo + 0.1, 5)
    expect_equal(.bnd$theta[["wt.cl"]], .lo)
  })
})
