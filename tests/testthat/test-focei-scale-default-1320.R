nmTest({
  # the FOCEi scale setup's count of thetas the linear rule scales
  .nLinear <- function(mod) {
    .ui <- rxode2::rxUiDecompress(rxode2::rxode2(mod))
    .th <- .ui$iniDf[!is.na(.ui$iniDf$ntheta), , drop = FALSE]
    .env <- new.env(parent = emptyenv())
    .env$lower <- .th$lower[order(.th$ntheta)]
    .foceiOptEnvSetupScaleC(.ui, .env)
    .env$nLinearTheta
  }

  test_that("the linear-theta count skips exp() thetas, residuals and zero starts", {
    .logOnly <- function() {
      ini({
        tcl <- 1.5
        eta.cl ~ 0.1
        add.sd <- 0.3
      })
      model({
        y <- exp(tcl + eta.cl)
        y ~ add(add.sd)
      })
    }
    .mixed <- function() {
      ini({
        tcl <- 1.5
        b <- 2.5
        bZero1320 <- 0
        eta.cl ~ 0.1
        add.sd <- 0.3
      })
      model({
        y <- exp(tcl + eta.cl) + b * x1 + bZero1320 * x2
        y ~ add(add.sd)
      })
    }
    expect_equal(.nLinear(.logOnly), 0L)
    # b is linear; a zero start keeps its zero-theta scale and is not counted
    nlmixr2global$nlmixr2EstEnv$zeroThetaScale <- data.frame(name = "bZero1320", est = 0, scale = 1)
    on.exit(nlmixr2global$nlmixr2EstEnv$zeroThetaScale <- NULL, add = TRUE)
    expect_equal(.nLinear(.mixed), 1L)
  })

  test_that("by default FOCEi adapts the scale only when a theta is linear", {
    skip_on_cran()
    rxode2::setRxThreads(1L)
    on.exit(rxode2::setRxThreads(NULL), add = TRUE)
    .ctl <- foceiControl(print = 0, covMethod = "", calcTables = FALSE)
    .theo <- function() {
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
    .f <- suppressMessages(suppressWarnings(nlmixr2(.theo, nlmixr2data::theo_sd, est = "focei", control = .ctl)))
    # no linear theta: no rescaled restart was even tried
    expect_null(.f$env$optReturn$nScaleAdapt)

    # the #1320 Emax model: two linear thetas, |init| scale and rescaled restarts
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
    .emax <- function() {
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
    .f <- suppressMessages(suppressWarnings(nlmixr2(.emax, .d, est = "focei", control = .ctl)))
    expect_equal(.f$scaleInfo$scaleC[1:2], c(9, 6))
    expect_false(is.null(.f$env$optReturn$nScaleAdapt))
    expect_lt(.f$objf, 165)
  })
})
