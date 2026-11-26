.catCovMod <- function() {
  ini({
    tcl <- log(10)
    tv <- log(50)
    tka <- log(1)
    clSex <- -0.1
    propSd <- c(0, 0.3, 1)
    eta.cl ~ 0.1
    eta.v ~ 0.1
    eta.ka ~ 0.1
  })
  model({
    cl <- exp(tcl + eta.cl) * (1 + clSex * SEX)
    v <- exp(tv + eta.v)
    ka <- exp(tka + eta.ka)
    linCmt() ~ prop(propSd)
  })
}

.catCovData <- function(n = 60, clSex = -0.25) {
  .sim <- .catCovMod |>
    ini(
      tcl = log(14.1),
      tv = log(62.5),
      tka = log(0.66),
      clSex = clSex,
      propSd = 0.2,
      eta.cl ~ 0.134,
      eta.v ~ 0.134,
      eta.ka ~ 0.16
    )
  .obsT <- c(0.5, 1, 2, 4, 8, 11.5, 60.5, 61, 62, 64, 68, 72)
  .ev <- do.call(
    rbind,
    lapply(seq_len(n), function(i) {
      rbind(
        data.frame(id = i, time = (0:5) * 12, amt = 450, evid = 1L, cmt = "depot"),
        data.frame(id = i, time = .obsT, amt = 0, evid = 0L, cmt = "central")
      )
    })
  )
  .ev$SEX <- rep(c(0, 1), length.out = n)[.ev$id]
  .s <- rxode2::rxWithSeed(1315, {
    rxode2::rxSolve(.sim, .ev, returnType = "data.frame", addDosing = TRUE, keep = "SEX")
  })
  data.frame(
    ID = .s$id,
    TIME = .s$time,
    AMT = ifelse(is.na(.s$amt), 0, .s$amt),
    EVID = .s$evid,
    CMT = ifelse(.s$evid == 0, 2L, 1L),
    DV = ifelse(.s$evid == 0, .s$sim, 0),
    SEX = .s$SEX
  )
}

test_that("a 0/1 multiplicative covariate is mu-referenced for saem only", {
  .d <- data.frame(ID = 1:4, TIME = 0, DV = 1, SEX = c(0, 1, 0, 1))
  .ctl <- saemControl()
  .mod <- rxode2::rxode2(.catCovMod)
  .r <- suppressWarnings(.preProcessSaemCatCov(.mod, "saem", .d, .ctl))
  expect_type(.r, "list")
  .ui <- .r$ui
  expect_true("rxBoundedTr.clSex" %in% .ui$muRefCovariateDataFrame$covariateParameter)
  expect_equal(.ui$iniDf$est[.ui$iniDf$name == "rxBoundedTr.clSex"], log1p(-0.1))
  .tr <- .ui$boundedTransforms
  expect_length(.tr, 1L)
  expect_identical(.tr[[1]]$type, "log1p")
  expect_identical(.tr[[1]]$lineOrig, "cl <- exp(tcl + eta.cl) * (1 + clSex * SEX)")
  expect_warning(
    .preProcessSaemCatCov(.mod, "saem", .d, .ctl),
    "'clSex' fit as log(1+clSex), mu-referenced",
    fixed = TRUE
  )
  # other methods, the control switch, a non-binary covariate
  expect_null(.preProcessSaemCatCov(.mod, "focei", .d, .ctl))
  expect_null(.preProcessSaemCatCov(.mod, "saem", .d, saemControl(muRefCovAlg = FALSE)))
  .d2 <- .d
  .d2$SEX <- c(0, 1, 2, 1)
  expect_null(.preProcessSaemCatCov(.mod, "saem", .d2, .ctl))
  # a theta used elsewhere cannot move into the mu-referenced phi
  .m2 <- .mod |> model(v <- exp(tv + eta.v) * (1 + clSex))
  expect_null(.preProcessSaemCatCov(.m2, "saem", .d, .ctl))
  # factor order and the covariate written first
  .m3 <- .mod |> model(cl <- (1 + SEX * clSex) * exp(tcl + eta.cl))
  .r3 <- suppressWarnings(.preProcessSaemCatCov(.m3, "saem", .d, .ctl))
  expect_true("rxBoundedTr.clSex" %in% .r3$ui$muRefCovariateDataFrame$covariateParameter)
  # runs after the other rewrites, just before the bounded transform
  .h <- preProcessHooks()
  expect_identical(tail(.h, 2L), c(".preProcessSaemCatCov", ".preProcessBoundedTransform"))
})

test_that("log1p back-transform and bounds", {
  .tr <- list(type = "log1p", lower = -0.99, upper = 0.5)
  expect_equal(.boundedTransformBack(.tr, log1p(-0.25)), -0.25)
  expect_equal(.boundedTransformDeriv(.tr, -0.25), 0.75)
  expect_equal(.boundedTransformBounds(.tr, -0.25), c(-0.99, 0.5))
  expect_equal(.boundedTransformBounds(.tr, 0.8), c(-0.99, 0.8))
  expect_equal(.boundedTransformBounds(list(type = "logit", lower = 0, upper = 1), 2), c(0, 1))
})

test_that("saem recovers the sign of a multiplicative 0/1 covariate (#1315)", {
  skip_on_cran()
  .d <- .catCovData()
  .ctl <- saemControl(print = 0, nBurn = 150, nEm = 150)
  .fit <- suppressMessages(nlmixr2(.catCovMod, .d, "saem", control = .ctl))
  expect_lt(.fit$theta[["clSex"]], -0.1)
  expect_gt(.fit$theta[["clSex"]], -0.45)
  expect_true(any(grepl("'clSex' fit as log(1+clSex)", .fit$runInfo, fixed = TRUE)))
  # the fit reports the model as written
  expect_false("rxBoundedTr.clSex" %in% .fit$iniDf$name)
  expect_true(any(vapply(
    .fit$ui$lstExpr,
    identical,
    logical(1),
    quote(cl <- exp(tcl + eta.cl) * (1 + clSex * SEX))
  )))
  expect_true("clSex" %in% names(.fit$parHist))
  expect_true("clSex" %in% rownames(.fit$cov))
  # delta-method SE from log(1 + clSex)
  expect_gt(.fit$parFixedDf["clSex", "SE"], 0)

  # bounds (and another bounded theta) do not stop the rewrite
  .mb <- .catCovMod |> ini(clSex = c(-0.99, -0.1, 5), tka = c(-5, 0, 5))
  .fb <- suppressMessages(suppressWarnings(
    nlmixr2(.mb, .d, "saem", control = saemControl(print = 0, nBurn = 50, nEm = 50))
  ))
  expect_lt(.fb$theta[["clSex"]], -0.1)
  expect_true(any(grepl("bounds of 'clSex' not enforced", .fb$runInfo, fixed = TRUE)))
  expect_equal(.fb$iniDf$lower[.fb$iniDf$name == "clSex"], -0.99)
})
