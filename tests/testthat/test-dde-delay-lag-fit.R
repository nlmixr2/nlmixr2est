# FOCEi recovers a delay() duration from a perturbed start (#1300).  Before the
# fix the lag's eta sensitivity was zero and lTR stalled ~22 OFV points short.

nmTest({
  .ddeLagMod <- function() {
    ini({
      lkel <- log(0.3)
      lv <- log(5)
      lkin <- log(2)
      lTR <- log(6)
      ls <- log(0.5)
      eta.kin ~ 0.05
      eta.TR ~ 0.05
      pkadd <- 0.1
      pdprop <- 0.05
    })
    model({
      kel <- exp(lkel)
      v <- exp(lv)
      kin <- exp(lkin + eta.kin)
      TR <- exp(lTR + eta.TR)
      s <- exp(ls)
      C <- central / v
      R(0) <- kin * TR
      d/dt(central) <- -kel * central
      d/dt(R) <- kin * s * C - kin * s * delay(central, TR) / v
      C ~ add(pkadd) | central
      R ~ prop(pdprop) | R
    })
  }

  .ddeLagData <- function() {
    .testSeed(7)
    .ev <- rxode2::et(amt = 50, cmt = "central") |>
      rxode2::et(c(0.5, 2, 4, 6, 8, 10, 12, 16, 24)) |>
      rxode2::et(id = 1:12)
    .sim <- rxode2::rxSolve(.ddeLagMod, .ev, returnType = "data.frame", addDosing = TRUE)
    .obs <- .sim$evid == 0
    .pd <- .sim[, c("id", "time", "evid", "amt")]
    .pd$cmt <- ifelse(.obs, "R", "central")
    .pd$dv <- ifelse(.obs, .sim$R * (1 + stats::rnorm(nrow(.sim), 0, 0.05)), NA_real_)
    .pk <- .sim[.obs, c("id", "time", "evid", "amt")]
    .pk$cmt <- "central"
    .pk$dv <- .sim$C[.obs] + stats::rnorm(sum(.obs), 0, 0.1)
    .d <- rbind(.pd, .pk)
    .d[order(.d$id, .d$time, -.d$evid), ]
  }

  test_that("the delay() duration's eta and theta are finite-differenced", {
    skip_on_cran()
    .ui <- rxode2::rxode2(.ddeLagMod)
    .m <- .ui$foceiModel
    expect_equal(.m$eventEta, c(0L, 1L))
    expect_equal(which(.m$eventTheta == 1L), 4L)
  })

  test_that("focei moves a delay() duration off its initial estimate", {
    skip_on_cran()
    .d <- .ddeLagData()
    .ctl <- foceiControl(print = 0L, covMethod = "")
    .truth <- .nlmixr(.ddeLagMod, .d, "focei", control = .ctl)
    .pert <- .nlmixr(.ddeLagMod |> rxode2::ini(lTR = log(10)), .d, "focei", control = .ctl)
    expect_lt(.pert$objf, .truth$objf + 2)
    expect_equal(unname(fixef(.pert)["lTR"]), log(6), tolerance = 0.15)
  })

  test_that("the analytic covariance declines a delay() duration parameter", {
    skip_on_cran()
    .ui <- rxode2::rxode2(.ddeLagMod)
    expect_null(.foceiAnalyticAugModelDirs(.ui, c("THETA_4_", "ETA_2_"), order = 1L))
    expect_match(.foceiAnalyticWhy$detail, "delay")
  })
})
