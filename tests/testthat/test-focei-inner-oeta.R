test_that("innerOeta is optimized on demand, once (#1321)", {
  skip_on_cran()
  model <- function() {
    ini({
      tka <- 0.45
      tcl <- 1
      tv <- 3.45
      prop.sd <- 0.1
      eta.ka ~ 0.6
      eta.cl ~ 0.3
    })
    model({
      ka <- exp(tka + eta.ka)
      cl <- exp(tcl + eta.cl)
      v <- exp(tv)
      d / dt(depot) <- -ka * depot
      d / dt(center) <- ka * depot - cl / v * center
      cp <- center / v + 0 * t
      cp ~ prop(prop.sd)
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
  expect_null(attributes(oeta))
  expect_identical(readRDS(cacheFile)$innerOetaOpt, oeta)
  # the raw innerOeta the bundle already carried is left alone
  expect_identical(suppressMessages(rxUiGet.foceiModel(ui))$innerOeta, raw)

  # a second request (and a reloaded bundle) reuses the stored text
  expect_identical(rxUiGet.foceiInnerOeta(ui), oeta)
  expect_equal(nOpt, 1L)

  # solvable as is, with the inner model's compartment order
  oetaMod <- rxode2::rxode2(oeta)
  expect_identical(rxode2::rxModelVars(oetaMod)$state, rxode2::rxModelVars(m$inner)$state)

  # the optimized and raw text solve to the same outputs
  wrap <- attr(raw, "nlmixr2estWrap")
  rawMod <- rxode2::rxode2(paste(c(wrap[1], raw, wrap[2]), collapse = "\n"))
  ev <- rxode2::et(amt = 100, cmt = "depot")
  ev <- rxode2::et(ev, c(0.5, 1, 2, 4, 8))
  p <- c(
    "THETA[1]" = 0.45,
    "THETA[2]" = 1,
    "THETA[3]" = 3.45,
    "THETA[4]" = 0.1,
    "ETA[1]" = 0.1,
    "ETA[2]" = -0.2
  )
  s1 <- as.data.frame(rxode2::rxSolve(rawMod, ev, p))
  s2 <- as.data.frame(rxode2::rxSolve(oetaMod, ev, p))
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
  expect_true(all(s2$rx_pred_ > 0))
})

test_that("foceiInnerOeta keeps the model's mtime() lines (#1321)", {
  skip_on_cran()
  model <- function() {
    ini({
      tke <- -1
      tt1 <- log(2)
      add.sd <- 0.1
      eta.ke ~ 0.1
    })
    model({
      ke <- exp(tke + eta.ke)
      mtime(t1) <- exp(tt1)
      cp <- 10 * exp(-ke * t)
      if (t >= t1) cp <- 5 * exp(-ke * t)
      cp ~ add(add.sd)
    })
  }
  ui <- list(rxode2::.copyUi(suppressMessages(nlmixr2(model))))
  cacheFile <- rxUiGet.foceiModelCache(ui)
  unlink(cacheFile)
  on.exit(unlink(cacheFile), add = TRUE)
  oeta <- suppressMessages(rxUiGet.foceiInnerOeta(ui))
  expect_match(oeta, "mtime(t1)", fixed = TRUE)
  oetaMod <- rxode2::rxode2(oeta)
  p <- c("THETA[1]" = -1, "THETA[2]" = log(2), "THETA[3]" = 0.1, "ETA[1]" = 0)
  s <- as.data.frame(rxode2::rxSolve(oetaMod, rxode2::et(c(1, 3)), p))
  s <- s[s$time %in% c(1, 3), ]
  expect_equal(s$rx_pred_, c(10 * exp(-exp(-1)), 5 * exp(-3 * exp(-1))), tolerance = 1e-6)
})
