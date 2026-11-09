test_that("innerOeta is optimized on demand, once (#1321)", {
  skip_on_cran()
  model <- function() {
    ini({
      tka <- 0.45
      tcl <- 1
      tv <- 3.45
      add.sd <- 0.7
      eta.ka ~ 0.6
      eta.cl ~ 0.3
    })
    model({
      ka <- exp(tka + eta.ka)
      cl <- exp(tcl + eta.cl)
      v <- exp(tv)
      cp <- 100 * ka / (v * (ka - cl / v)) * (exp(-cl / v * t) - exp(-ka * t))
      cp ~ add(add.sd)
    })
  }
  ui <- list(rxode2::.copyUi(suppressMessages(nlmixr2(model))))
  cacheFile <- rxUiGet.foceiModelCache(ui)
  unlink(cacheFile)
  on.exit(unlink(cacheFile), add = TRUE)

  # count top-level passes only; rxOptExpr() re-enters itself per chunk
  nOpt <- 0L
  inOpt <- FALSE
  realOpt <- rxode2::rxOptExpr
  local_mocked_bindings(
    rxOptExpr = function(x, msg = "model", ...) {
      if (!inOpt && any(grepl("rx__ETA1=", x, fixed = TRUE))) {
        nOpt <<- nOpt + 1L
        inOpt <<- TRUE
        on.exit(inOpt <<- FALSE)
      }
      realOpt(x, msg, ...)
    },
    .package = "rxode2"
  )

  m <- suppressMessages(rxUiGet.foceiModel(ui))
  expect_equal(nOpt, 0L)
  raw <- m$innerOeta
  expect_true(is.character(raw))
  expect_false(is.null(attr(raw, "nlmixr2estOptimize")))

  oeta <- rxUiGet.foceiInnerOeta(ui)
  expect_equal(nOpt, 1L)
  expect_null(attr(oeta, "nlmixr2estOptimize"))
  expect_null(attr(readRDS(cacheFile)$innerOeta, "nlmixr2estOptimize"))

  # a second request (and a reloaded bundle) reuses the stored text
  expect_identical(rxUiGet.foceiInnerOeta(ui), oeta)
  expect_equal(nOpt, 1L)

  # the optimized and raw text solve to the same outputs
  d <- data.frame(id = 1, time = c(0.5, 1, 2, 4, 8), dv = 1)
  p <- c(
    "THETA[1]" = 0.45,
    "THETA[2]" = 1,
    "THETA[3]" = 3.45,
    "THETA[4]" = 0.7,
    "ETA[1]" = 0.1,
    "ETA[2]" = -0.2
  )
  rawTxt <- .addMtimeLines(raw, list(..mtime = attr(raw, "nlmixr2estMtime")))
  s1 <- as.data.frame(rxode2::rxSolve(rxode2::rxode2(rawTxt), d, p))
  s2 <- as.data.frame(rxode2::rxSolve(rxode2::rxode2(oeta), d, p))
  cols <- c(
    "rx_pred_",
    "rx_r_",
    "rx__ETA1",
    "rx__ETA2",
    "rx__sens_rx_pred__BY_ETA_1___",
    "rx__sens_rx_pred__BY_ETA_2___"
  )
  expect_equal(s1[, cols], s2[, cols], tolerance = 1e-8)
  expect_equal(s1$rx__ETA1, rep(0.1, 5))
})
