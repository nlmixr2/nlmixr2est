nmTest({
  .nlmeResidDat <- function(sdFun, n = 40) {
    rxode2::rxWithSeed(1464, {
      .t <- c(0.5, 1, 2, 4, 6, 8, 12, 16, 24)
      .d <- expand.grid(TIME = .t, ID = seq_len(n))
      .wt <- round(runif(n, 50, 110))
      .d$WT <- .wt[.d$ID]
      .ke <- 0.2 * exp(rnorm(n, 0, 0.3))
      .ipre <- 10 * exp(-.ke[.d$ID] * .d$TIME)
      .d$DV <- .ipre + rnorm(nrow(.d)) * sdFun(.ipre, .d$TIME, .d$WT)
      .d[, c("ID", "TIME", "WT", "DV")]
    })
  }

  # every ini() parameter the residual may use; the unused ones are dropped
  .nlmeResidIni <- list(
    add.sd = quote(add.sd <- 0.5),
    prop.sd = quote(prop.sd <- 0.2),
    pw = quote(pw <- 1),
    k = quote(k <- -0.1),
    eta.sd = quote(eta.sd ~ 0.1)
  )

  .nlmeResidModel <- function(...) {
    .err <- list(...)
    .vars <- unique(unlist(lapply(.err, all.vars)))
    .ini <- unname(.nlmeResidIni[names(.nlmeResidIni) %in% .vars])
    .f <- eval(bquote(
      function() {
        ini({
        tke <- log(0.2)
        eta.ke ~ 0.09
        ..(.ini)
      })
        model({
        ke <- exp(tke + eta.ke)
        ipre <- 10 * exp(-ke * t)
        ..(.err)
      })
      },
      splice = TRUE
    ))
    # rxode2 reads ini() labels from the source, which still holds `..()`
    attr(.f, "srcref") <- NULL
    suppressMessages(rxode2::rxode2(.f))
  }

  test_that("time- and covariate-scaled residual SDs fit under nlme (issue 1464)", {
    skip_on_cran()
    .cases <- list(
      addTime = list(
        quote(ipre ~ add(add.sd * exp(-time / 8))),
        function(f, t, wt) 0.3 * exp(-t / 8),
        "add.sd",
        0.3
      ),
      addWt = list(
        quote(ipre ~ add(add.sd * WT / 70)),
        function(f, t, wt) 0.3 * wt / 70,
        "add.sd",
        0.3
      ),
      propTime = list(
        quote(ipre ~ prop(prop.sd * (1 + time / 12))),
        function(f, t, wt) 0.1 * (1 + t / 12) * f,
        "prop.sd",
        0.1
      )
    )
    for (.n in names(.cases)) {
      .c <- .cases[[.n]]
      .dat <- .nlmeResidDat(.c[[2]])
      fit <- .nlmixr(.nlmeResidModel(.c[[1]]), .dat, "nlme", control = nlmeControl(verbose = FALSE))
      .th <- fit$theta
      expect_equal(.th[[.c[[3]]]], .c[[4]], tolerance = 0.15, info = .n)
      expect_equal(exp(.th[["tke"]]), 0.2, tolerance = 0.15, info = .n)
      # the theta reproduces nlme's own residual sd of every observation
      .nl <- fit$nlme
      .sdNlme <- as.numeric(residuals(.nl, type = "response") / residuals(.nl, type = "pearson"))
      .f <- as.numeric(fitted(.nl))
      .sd <- .c[[2]](.f, .dat$TIME, .dat$WT) / .c[[4]] * .th[[.c[[3]]]]
      expect_equal(.sd, .sdNlme, tolerance = 1e-6, info = .n)
      expect_false(any(c("add.sd", "prop.sd") %in% names(nlme::fixef(.nl))), info = .n)
    }
  })

  test_that("a time-scaled pow() residual SD fits under nlme (issue 1464)", {
    skip_on_cran()
    .sdFun <- function(f, t, wt) 0.1 * (1 + t / 12) * f^0.8
    .dat <- .nlmeResidDat(.sdFun)
    .m <- .nlmeResidModel(quote(ipre ~ pow(prop.sd * (1 + time / 12), pw)))
    fit <- .nlmixr(.m, .dat, "nlme", control = nlmeControl(verbose = FALSE))
    .th <- fit$theta
    expect_equal(.th[["prop.sd"]], 0.1, tolerance = 0.2)
    expect_equal(.th[["pw"]], 0.8, tolerance = 0.2)
    .nl <- fit$nlme
    .sdNlme <- as.numeric(residuals(.nl, type = "response") / residuals(.nl, type = "pearson"))
    .f <- as.numeric(fitted(.nl))
    .sd <- .th[["prop.sd"]] * (1 + .dat$TIME / 12) * abs(.f)^.th[["pw"]]
    expect_equal(.sd, .sdNlme, tolerance = 1e-6)
  })

  test_that("a fixed pow() exponent stays fixed under nlme (issue 1464)", {
    skip_on_cran()
    .dat <- .nlmeResidDat(function(f, t, wt) 0.1 * (1 + t / 12) * f^0.8)
    .m <- .nlmeResidModel(quote(ipre ~ pow(prop.sd * (1 + time / 12), pw)))
    .m <- suppressMessages(rxode2::ini(.m, pw = fix(0.8)))
    .info <- .nlmeResidInfo(.m)
    expect_equal(.info$pow, list(name = "pw", ini = 0.8, fix = TRUE))
    fit <- .nlmixr(.m, .dat, "nlme", control = nlmeControl(verbose = FALSE))
    expect_equal(fit$theta[["pw"]], 0.8)
    expect_equal(fit$theta[["prop.sd"]], 0.1, tolerance = 0.15)
  })

  test_that("an intermediate variable of the residual SD is resolved (issue 1464)", {
    .m <- .nlmeResidModel(quote(wn <- WT / 70), quote(ipre ~ add(add.sd * wn)))
    .info <- .nlmeResidInfo(.m)
    expect_equal(.info$theta, "add.sd")
    expect_equal(.info$covs, "WT")
    expect_false("add.sd" %in% .nlmeParamsToEstimate(list(.m)))
    expect_false("add.sd" %in% names(.m$nlmeStart))
  })

  test_that("unsupported residual errors are refused with a specific message (issue 1464)", {
    .refuse <- function(err, msg) {
      expect_error(.nlmeResidInfo(.nlmeResidModel(err)), msg, info = deparse1(err))
    }
    .refuse(quote(ipre ~ add(add.sd * exp(eta.sd))), "eta on the residual error \\('eta.sd'\\)")
    .refuse(quote(ipre ~ add(add.sd * ke)), "model-calculated 'ke'")
    .refuse(quote(ipre ~ add(add.sd * WT / 70) + prop(prop.sd)), "combined residual error")
    .refuse(quote(ipre ~ add(add.sd * exp(k * time))), "more than one estimated parameter")
    .refuse(quote(ipre ~ add(add.sd * exp(tke))), "sharing 'tke' with the structural model")
    .refuse(quote(ipre ~ add(add.sd * WT / 70)), NA)
    .m <- .nlmeResidModel(quote(ipre ~ pow(prop.sd * (1 + time / 12), pw)))
    .m <- suppressMessages(rxode2::ini(.m, prop.sd = fix(0.1)))
    expect_error(.nlmeResidInfo(.m), "without an estimated scale parameter")
  })

  test_that("a residual SD that is not proportional to its theta is refused (issue 1464)", {
    .m <- .nlmeResidModel(quote(ipre ~ add(add.sd^2 * WT / 70)))
    .info <- .nlmeResidInfo(.m)
    .dat <- data.frame(TIME = 1:3, WT = c(60, 70, 80))
    expect_error(.nlmeResidVarCov(.info, .dat), "must be 'add.sd' times")
    .m <- .nlmeResidModel(quote(ipre ~ add(add.sd * WT / 70)))
    expect_equal(.nlmeResidVarCov(.nlmeResidInfo(.m), .dat), (.dat$WT / 70)^2)
  })

  test_that("the nlme refusal message reaches nlmixr2 (issue 1464)", {
    .m <- .nlmeResidModel(quote(ipre ~ add(add.sd * exp(eta.sd))))
    expect_error(
      .nlmixr(.m, .nlmeResidDat(function(f, t, wt) 0.3), "nlme", control = nlmeControl(verbose = FALSE)),
      "nlme cannot estimate an eta on the residual error"
    )
  })
})
