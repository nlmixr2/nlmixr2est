nmTest({
  # one linear theta per initial value, across both edges of scaleCband
  .inits <- c(0.05, 0.1, 0.5, 5, 10, 20)
  .sweep <- eval(parse(
    text = paste0(
      "function() {\n ini({\n",
      paste0("  a", seq_along(.inits), " <- ", .inits, collapse = "\n"),
      "\n  eta.k ~ 0.1\n  sdw <- 1\n })\n model({\n  y <- eta.k + ",
      paste0("a", seq_along(.inits), " * x", seq_along(.inits), collapse = " + "),
      "\n  y ~ add(sdw)\n })\n}\n"
    )
  ))
  .scaleLinear <- function(mod, ...) {
    .ui <- rxode2::rxUiDecompress(rxode2::rxode2(mod))
    .ctl <- foceiControl(...)
    for (.n in c("scaleCband", "scaleLinear", "scaleC")) {
      rxode2::rxAssignControlValue(.ui, .n, .ctl[[.n]])
    }
    .ui$scaleCtheta
  }

  test_that("scaleLinear validates and round-trips", {
    expect_equal(foceiControl()$scaleLinear, "inverse")
    expect_equal(foceiControl(scaleLinear = "abs")$scaleLinear, "abs")
    expect_error(foceiControl(scaleLinear = "log"))
  })

  test_that("scaleLinear = 'inverse' keeps the banded 1/|init| default", {
    .sc <- .scaleLinear(.sweep)[seq_along(.inits)]
    expect_equal(.sc, c(0.05, 10, 2, 0.2, 0.1, 20))
  })

  test_that("scaleLinear = 'abs' scales every linear theta by |init|", {
    .sc <- .scaleLinear(.sweep, scaleLinear = "abs")
    expect_equal(.sc[seq_along(.inits)], .inits)
    # the residual error keeps its own 0.5 * |init|
    expect_equal(.sc[length(.sc)], 0.5)
  })

  test_that("scaleLinear = 'abs' takes |init| of a negative init, 1 at 0, and keeps a given scaleC", {
    .mod <- function() {
      ini({
        a <- -50
        bZero1320 <- 0
        c <- 3
        eta.k ~ 0.1
        add.sd <- 0.3
      })
      model({
        y <- eta.k + a * x1 + bZero1320 * x2 + c * x3
        y ~ add(add.sd)
      })
    }
    expect_equal(.scaleLinear(.mod, scaleLinear = "abs"), c(50, 1, 3, 0.15))
    expect_equal(.scaleLinear(.mod, scaleLinear = "abs", scaleC = 9), c(9, 1, 3, 0.15))
  })

  test_that("scaleLinear = 'abs' leaves mu-referenced thetas alone", {
    .mod <- function() {
      ini({
        tcl <- 1.5
        b <- 2.5
        eta.cl ~ 0.1
        add.sd <- 0.3
      })
      model({
        cl <- exp(tcl + eta.cl)
        y <- cl + b * x
        y ~ add(add.sd)
      })
    }
    expect_equal(.scaleLinear(.mod), c(1, 0.4, 0.15))
    expect_equal(.scaleLinear(.mod, scaleLinear = "abs"), c(1, 2.5, 0.15))
  })

  test_that("scaleLinear = 'abs' moves the stalled linear thetas of #1320", {
    skip_on_cran()
    rxode2::setRxThreads(1L)
    on.exit(rxode2::setRxThreads(NULL), add = TRUE)
    # arm-level Emax meta-analysis from #1320: tve0/tvemax start at 9 and 6,
    # inside scaleCband, so the default first steps are 0.25% and 0.6%
    set.seed(20261108)
    .arms <- expand.grid(ARM = 1:3, ID = 1:8)
    .arms$DOSE <- c(0, 50, 200)[.arms$ARM]
    .arms$N <- sample(40:200, nrow(.arms), replace = TRUE)
    .arms$OCC <- .arms$ARM
    .d <- merge(.arms, data.frame(TIME = c(0, 2, 4, 8, 12)))
    .d <- .d[order(.d$ID, .d$OCC, .d$TIME), ]
    .e0 <- 10 + stats::rnorm(8, 0, 1)[.d$ID]
    .emaxArm <- 8 * exp(stats::rnorm(nrow(.arms), 0, 0.1))[match(paste(.d$ID, .d$OCC), paste(.arms$ID, .arms$OCC))]
    .mean <- .e0 - .emaxArm * .d$DOSE / (40 + .d$DOSE) * .d$TIME / (3 + .d$TIME)
    .d$W <- .d$N / 9
    .d$DV <- (.mean + stats::rnorm(nrow(.d), 0, 1 / sqrt(.d$W))) * sqrt(.d$W)
    .mod <- function() {
      ini({
        tve0 <- 9
        tvemax <- 6
        led50 <- log(30)
        lt50 <- log(2)
        eta.e0 ~ 1
        iov.emax ~ 0.05 | OCC
        sdw <- 1
      })
      model({
        E0 <- tve0 + eta.e0
        EMAX <- tvemax * exp(iov.emax)
        ONSET <- 0
        if (TIME > 0) ONSET <- TIME / (exp(lt50) + TIME)
        EFF <- EMAX * DOSE / (exp(led50) + DOSE) * ONSET
        pred <- (E0 - EFF) * sqrt(W)
        pred ~ add(sdw)
      })
    }
    .ctl <- function(...) foceiControl(print = 0, covMethod = "", calcTables = FALSE, sigdig = 3, ...)
    .abs <- suppressMessages(suppressWarnings(
      nlmixr2(.mod, .d, est = "focei", control = .ctl(scaleLinear = "abs"))
    ))
    # the mechanism: the optimizer saw |init| for the linear thetas
    expect_equal(.abs$scaleInfo$scaleC[1:2], c(9, 6))
    # and reached the optimum FeRx and the band-free scale agree on (161.24)
    expect_lt(.abs$objf, 165)
    expect_gt(.abs$theta[["tvemax"]], 7.5)
  })
})
