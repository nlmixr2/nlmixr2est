# Event-modifier (f()/alag()) eta sensitivities -- nlmixr2est#1016.
#
# The reported symptom was a fit that left every f()/alag() eta at ~1e-9 while
# their omegas stayed finite, deterministically on the SECOND and later fits of
# a model in one R session.  The cause was the focei model bundle: it is cached
# as model TEXT and rehydrated with rxode2(), and eventSens is not recoverable
# from a compiled model, so the first fit built the jump models and every later
# one rebuilt them in "fd" mode.  Nothing errored -- the dosing-parameter
# sensitivity was simply zero, so those etas had no reason to move.
#
# What these tests pin, from the bottom up:
#   1. the cached-bundle round trip keeps the mode AND the shape still installs;
#   2. the inner model's d(pred)/d(eta) for a dosing eta matches central
#      differences (so the jumps are not merely present but right);
#   3. a repeat fit in one session gives the same etas as the first, and the
#      dosing etas are not collapsed.
# (1) and (3) are what regress if the mode is ever dropped again.

nmTest({

  .evEtaMod <- function() {
    ini({
      tka <- log(1.2)
      tcl <- log(2)
      tv <- log(20)
      tf <- -0.5
      tlag <- log(0.5)
      eta.cl ~ 0.1
      eta.f ~ 0.1
      eta.lag ~ 0.1
      add.sd <- 0.5
    })
    model({
      ka <- exp(tka)
      cl <- exp(tcl) * exp(eta.cl)
      v <- exp(tv)
      f(depot) <- expit(tf + eta.f)
      alag(depot) <- exp(tlag + eta.lag)
      d/dt(depot) <- -ka * depot
      d/dt(central) <- ka * depot - (cl / v) * central
      cp <- central / v
      cp ~ add(add.sd)
    })
  }

  # non-uniform dose/observation spacing (uniform spacing has hidden pairing
  # bugs in this area before)
  .evEtaDat <- function(nid = 4L) {
    do.call(rbind, lapply(seq_len(nid), function(i) {
      tim <- c(0, 3, 7, 15, 24, 30, 41, 50) + (i - 1) * 0.5
      data.frame(
        id = i, time = tim, amt = c(100, 0, 100, 0, 100, 0, 100, 0),
        evid = c(1, 0, 1, 0, 1, 0, 1, 0), cmt = 1,
        dv = c(0, 1.7, 2.4, 1.1, 2.9, 2.2, 3.1, 1.8) + 0.1 * i
      )
    }))
  }

  .evEtaBundle <- function(eventSens) {
    .ui <- rxode2::.copyUi(suppressMessages(nlmixr2est::nlmixr2(.evEtaMod)))
    assign("control", foceiControl(eventSens = eventSens), envir = .ui)
    suppressMessages(suppressWarnings(.ui$foceiModel))
  }

  test_that("a dosing eta is flagged in eventEtaAll in both eventSens modes", {
    skip_on_cran()
    # eta.cl is not a dosing parameter; eta.f and eta.lag are
    for (.es in c("jump", "fd")) {
      .b <- .evEtaBundle(.es)
      expect_equal(as.integer(.b$eventEtaAll), c(0L, 1L, 1L))
    }
    # the finite-difference switches C++ reads stay off under "jump" (the
    # analytic jump sensitivity is computed, so the FD fallback must not run)
    expect_equal(as.integer(.evEtaBundle("jump")$eventEta), c(0L, 0L, 0L))
    expect_equal(as.integer(.evEtaBundle("fd")$eventEta), c(0L, 1L, 1L))
  })

  test_that("the cached bundle round trip keeps the jump shape (#1016)", {
    skip_on_cran()
    .inner <- .evEtaBundle("jump")$inner
    expect_equal(attr(.inner, "nlmixr2estEventSens"), "jump")
    # deflate/inflate is what every fit after the first one does
    .rt <- .foceiModelCacheInflate(.foceiModelCacheDeflate(.inner))
    expect_equal(attr(.rt, "nlmixr2estEventSens"), "jump")
    # and it round trips again -- an inflated bundle must not deflate to NULL
    .rt2 <- .foceiModelCacheInflate(.foceiModelCacheDeflate(.rt))
    expect_equal(attr(.rt2, "nlmixr2estEventSens"), "jump")
    # the point of the mode: the event-sensitivity shape still installs.  This
    # is the exact call .foceiFitInternal makes before the C++ solve loop, and
    # it returned FALSE for a rehydrated bundle before the fix.
    on.exit(rxode2::rxEventSensDeactivate(), add = TRUE)
    expect_true(isTRUE(rxode2::rxEventSensLoadModel(.rt2)))
  })

  test_that("a jump fit that lost its shape says so (#1016)", {
    skip_on_cran()
    .b <- .evEtaBundle("jump")
    # the shape installed -- nothing to say
    expect_false(.foceiEventSensWarn(TRUE, .b))
    expect_silent(.foceiEventSensWarn(TRUE, .b))
    # an "fd" fit never installs a shape; its dosing etas go through the C++
    # finite-difference fallback, so a missing shape there is not a defect
    expect_false(.foceiEventSensWarn(FALSE, .b, "fd"))
    expect_silent(.foceiEventSensWarn(FALSE, .b, "fd"))
    # it did not, and this model's etas need it: the fit would otherwise leave
    # eta.f/eta.lag where they started with no diagnostic at all
    expect_warning(
      expect_true(.foceiEventSensWarn(FALSE, .b)),
      "event sensitivities not loaded"
    )
    # the warning reaches the user through $runInfo, which renders one bullet
    # per line
    expect_lt(
      nchar(tryCatch(.foceiEventSensWarn(FALSE, .b),
        warning = function(w) conditionMessage(w)
      )), 75
    )
    # a model with no dosing eta does not depend on the jumps, so it stays quiet
    .b$eventEtaAll <- c(0L, 0L, 0L)
    expect_false(.foceiEventSensWarn(FALSE, .b))
    expect_silent(.foceiEventSensWarn(FALSE, .b))
    # unknown is NOT quiet: the scan that fills eventEtaAll only runs on a model
    # that doses through f()/alag()/rate()/dur(), so a failed scan cannot be
    # read as "no dosing etas" -- that would silence the tripwire on exactly the
    # models it exists for
    .b$eventEtaAll <- c(NA_integer_, NA_integer_, NA_integer_)
    expect_warning(expect_true(.foceiEventSensWarn(FALSE, .b)))
    # a bundle from before the field existed establishes nothing either way
    .b$eventEtaAll <- NULL
    expect_false(.foceiEventSensWarn(FALSE, .b))
    expect_silent(.foceiEventSensWarn(FALSE, .b))
  })

  test_that("dosing-eta sensitivities match central differences (#1016)", {
    skip_on_cran()
    .m <- .evEtaBundle("jump")$inner
    .ev <- .evEtaDat(1L)
    .p <- c(
      `THETA[1]` = log(1.2), `THETA[2]` = log(2), `THETA[3]` = log(20),
      `THETA[4]` = -0.5, `THETA[5]` = log(0.5), `THETA[6]` = 0.5,
      `ETA[1]` = 0.3, `ETA[2]` = -0.2, `ETA[3]` = 0.1
    )
    .slv <- function(q) {
      suppressWarnings(rxode2::rxSolve(.m,
        params = q, events = .ev, returnType = "data.frame",
        atol = 1e-12, rtol = 1e-12
      ))
    }
    .r0 <- .slv(.p)
    .h <- 1e-4
    for (.k in 2:3) { # ETA[2] = f(), ETA[3] = alag()
      .en <- paste0("ETA[", .k, "]")
      .a <- .p
      .a[.en] <- .a[.en] + .h
      .b <- .p
      .b[.en] <- .b[.en] - .h
      .fd <- (.slv(.a)$rx_pred_ - .slv(.b)$rx_pred_) / (2 * .h)
      .an <- .r0[[paste0("rx__sens_rx_pred__BY_ETA_", .k, "___")]]
      # not merely nonzero: right.  Without the jumps .an is exactly 0.
      expect_gt(max(abs(.an)), 1e-3)
      expect_lt(max(abs(.an - .fd)) / max(abs(.fd)), 1e-4)
    }
  })

  test_that("a repeat fit keeps its f()/alag() etas (#1016)", {
    skip_on_cran()
    .dat <- .evEtaDat()
    .ctl <- foceiControl(
      print = 0, maxOuterIterations = 0L, covMethod = "",
      calcTables = FALSE, sigdig = 6, etaNudge = 0, etaNudge2 = 0
    )
    .f1 <- .nlmixr(.evEtaMod, .dat, est = "focei", control = .ctl)
    # the second fit of a model in one session is the one that read the cached
    # bundle -- it lost the jump mode, and with it every dosing-eta gradient
    .f2 <- .nlmixr(.evEtaMod, .dat, est = "focei", control = .ctl)
    expect_equal(.f2$objective, .f1$objective, tolerance = 1e-6)
    expect_equal(as.matrix(.f2$eta[, -1]), as.matrix(.f1$eta[, -1]),
      tolerance = 1e-5
    )
    # and they did not collapse: with the jumps gone both fits sat at ~1e-9
    expect_gt(max(abs(.f1$eta$eta.f)), 1e-3)
    expect_gt(max(abs(.f2$eta$eta.f)), 1e-3)
  })

  test_that("a fit whose cached bundle lost the mode says so (#1016)", {
    skip_on_cran()
    .dat <- .evEtaDat()
    .ctl <- foceiControl(
      print = 0, maxOuterIterations = 0L, covMethod = "",
      calcTables = FALSE, sigdig = 6, etaNudge = 0, etaNudge2 = 0
    )
    .msg <- "event sensitivities not loaded"
    # healthy first: the assertion below has to be able to come out either way
    .ok <- .nlmixr(.evEtaMod, .dat, est = "focei", control = .ctl)
    expect_false(any(grepl(.msg, .ok$runInfo)))

    # now reproduce #1016 itself end to end.  The bundle is cached as model
    # TEXT plus the mode that built it; before the fix the mode was not stored,
    # so every fit after the first rehydrated the inner model in "fd" mode.
    # Poisoning the stored mode is exactly that state.
    .ui <- rxode2::.copyUi(suppressMessages(nlmixr2est::nlmixr2(.evEtaMod)))
    assign("control", .ctl, envir = .ui)
    .cacheFile <- .ui$foceiModelCache
    on.exit(unlink(.cacheFile), add = TRUE)
    expect_true(file.exists(.cacheFile))
    .store <- readRDS(.cacheFile)
    expect_equal(.store$inner$eventSens, "jump")
    .store$inner$eventSens <- "fd"
    saveRDS(.store, .cacheFile)

    # the fit still runs and still returns -- that is the whole problem -- so
    # the only thing that tells anyone is the note on $runInfo
    .bad <- .nlmixr(.evEtaMod, .dat, est = "focei", control = .ctl)
    expect_true(any(grepl(.msg, .bad$runInfo)))
    # and it is the #1016 answer: the dosing etas did not move (reported as
    # ~9e-09 on the issue's own model; the healthy fit below is 100x larger)
    expect_lt(max(abs(.bad$eta$eta.f)), 1e-5)
    expect_gt(max(abs(.ok$eta$eta.f)), 1e-3)
  })

  # dosing built from THETAs alone: the tripwire is eta-scoped on purpose, and
  # this is why.  A missing inner shape cannot change such a fit -- the inner
  # model carries eta sensitivities only, a theta gradient comes from the outer
  # re-solve, and the bundle's eventTheta is read nowhere in src/ -- so warning
  # on it would be a false alarm on every model with a modeled lag and no eta
  # on it.
  .evThetaMod <- function() {
    ini({
      tka <- log(1.2)
      tcl <- log(2)
      tv <- log(20)
      tf <- -0.5
      tlag <- log(0.5)
      eta.cl ~ 0.1
      add.sd <- 0.5
    })
    model({
      ka <- exp(tka)
      cl <- exp(tcl) * exp(eta.cl)
      v <- exp(tv)
      f(depot) <- expit(tf)
      alag(depot) <- exp(tlag)
      d/dt(depot) <- -ka * depot
      d/dt(central) <- ka * depot - (cl / v) * central
      cp <- central / v
      cp ~ add(add.sd)
    })
  }

  test_that("dosing on thetas alone is untouched by a lost shape (#1016)", {
    skip_on_cran()
    .dat <- .evEtaDat()
    .ctl <- foceiControl(
      print = 0, maxOuterIterations = 20L, covMethod = "",
      calcTables = FALSE, sigdig = 5, etaNudge = 0, etaNudge2 = 0
    )
    .ui <- rxode2::.copyUi(suppressMessages(nlmixr2est::nlmixr2(.evThetaMod)))
    assign("control", .ctl, envir = .ui)
    .cacheFile <- .ui$foceiModelCache
    on.exit(unlink(.cacheFile), add = TRUE)
    .ok <- .nlmixr(.evThetaMod, .dat, est = "focei", control = .ctl)
    # no eta enters a dosing expression, so there is nothing for the jumps to
    # carry into the inner problem
    .store <- readRDS(.cacheFile)
    expect_equal(as.integer(.store$eventEtaAll), 0L)

    .store$inner$eventSens <- "fd"
    saveRDS(.store, .cacheFile)
    .bad <- .nlmixr(.evThetaMod, .dat, est = "focei", control = .ctl)
    # identical fit, to the last digit -- measured max abs theta difference 0
    expect_equal(.bad$objective, .ok$objective, tolerance = 1e-10)
    expect_equal(.bad$theta, .ok$theta, tolerance = 1e-10)
    # so it must stay quiet
    expect_false(any(grepl("event sensitivities not loaded", .bad$runInfo)))
  })

})
