nmTest({
  # shi21Forward()/shi21Central() (src/shi21.cpp): the Shi et al. (2021) step
  # search when one of its outer probes is not finite.

  # eventSens = "fd" finite-differences fdepot (a dosing parameter) for every
  # subject; the log-likelihood is -Inf once fdepot passes `edge` (it starts at
  # 0.9) by more than about 710 / k, so a probe past it is not finite.  (rxode2
  # reads a model from its source when it has one, which still says .(edge).)
  .shiMod <- function(edge, k = 1e9) {
    utils::removeSource(eval(bquote(function() {
      ini({
        fdepot <- 0.9
        tka <- 0.45
        tcl <- 1
        tv <- 3.45
        add.sd <- 0.7
      })
      model({
        ka <- exp(tka)
        cl <- exp(tcl)
        v <- exp(tv)
        d / dt(depot) <- -ka * depot
        f(depot) <- fdepot
        d / dt(center) <- ka * depot - cl / v * center
        cp <- center / v
        ll(err) ~ -0.5 * ((DV - cp) / add.sd)^2 - log(add.sd) - exp(.(k) * (fdepot - .(edge)))
      })
    })))
  }

  # The fdepot gradient of the first evaluation, whose one-iteration step search
  # (shi21maxFD = 1) gives it; the difference it should equal at the search's
  # starting step h0 (on the natural scale); and a central difference at 1e-6
  .shiFirstGrad <- function(edge, eventType) {
    .ctl <- nlmControl(print = 0L, eventType = eventType, shi21maxFD = 1L, eventSens = "fd")
    .x <- suppressMessages(nlmObjectiveSetup(
      .shiMod(edge),
      nlmixr2data::theo_sd,
      control = .ctl,
      gradient = TRUE,
      scale = "natural"
    ))
    on.exit(.nlmFreeEnv())
    .ef <- .Machine$double.eps^(1 / 3)
    .e <- c(1, 0, 0, 0, 0)
    if (eventType == "central") {
      .h0 <- (3 * .ef)^(1 / 3)
      .diff <- (nlmSolveR(.x + .h0 * .e) - nlmSolveR(.x - .h0 * .e)) / (2 * .h0)
      .outer <- nlmSolveR(.x + 3 * .h0 * .e)
    } else {
      .h0 <- 2 / sqrt(3) * sqrt(.ef)
      .diff <- (nlmSolveR(.x + .h0 * .e) - nlmSolveR(.x)) / .h0
      .outer <- nlmSolveR(.x + 4 * .h0 * .e)
    }
    .fine <- (nlmSolveR(.x + 1e-6 * .e) - nlmSolveR(.x - 1e-6 * .e)) / 2e-6
    c(grad = nlmLikEvalC_(.x)$grad[1], diff = .diff, fine = .fine, outer = .outer)
  }

  test_that("a step search whose outer probe is not finite keeps the difference at its step (issue 1140)", {
    skip_on_cran()
    # central: f(x +- h0) are finite and f(x + 3 h0) is not, so the gradient
    # is the central difference at h0 (the search then shrinks the step)
    .r <- .shiFirstGrad(0.9 + 0.05, "central")
    expect_false(is.finite(.r[["outer"]]))
    expect_equal(.r[["grad"]], .r[["diff"]], tolerance = 1e-6)
    # forward: f(x + h0) is finite and f(x + 4 h0) is not, so the gradient is
    # the forward difference at h0, within 1.5% of the derivative
    .r <- .shiFirstGrad(0.9 + 0.006, "forward")
    expect_false(is.finite(.r[["outer"]]))
    expect_equal(.r[["grad"]], .r[["fine"]], tolerance = 0.02)
    # f(x) is the pred model's too, the model f(x + h0) is solved with
    expect_equal(.r[["grad"]], .r[["diff"]], tolerance = 1e-6)
  })

  test_that("the step shi21Central() returns is the one its gradient differences (issue 1140)", {
    # x^2 at 1, finite only on [1, 1.01]: the first probe (1 + h0, h0 = 0.03)
    # fails, the search shrinks to h0 / 6 and takes a forward difference there
    # (the backward probe fails), and every later backward probe fails too, so
    # that is the difference and the step it returns
    .f <- function(x) if (x < 1 || x > 1.01) NA_real_ else x^2
    .r <- shi21CentralWrap(.f, 1, 1, 1L, 0.03^3 / 3)
    expect_equal(.r$h, 0.005, tolerance = 1e-12)
    # the forward difference of x^2 at 1 is 2 + h
    expect_equal(drop(.r$gr), 2 + .r$h, tolerance = 1e-12)
  })

  test_that("a theta whose plus-side probes never solve takes the backward difference", {
    skip_on_cran()
    withr::defer(.nlmFreeEnv())
    # the log-likelihood is -Inf at any fdepot above 0.9, an edge sharp enough
    # that every probe the search tries past it fails, so the fdepot gradient is
    # the backward difference at the search's starting step h0
    .mod <- .shiMod(0.9, 1e300)
    .ef <- .Machine$double.eps^(1 / 3)
    .e <- c(1, 0, 0, 0, 0)
    for (.et in c("central", "forward")) {
      .ctl <- nlmControl(print = 0L, eventType = .et, eventSens = "fd")
      .x <- suppressMessages(nlmObjectiveSetup(
        .mod,
        nlmixr2data::theo_sd,
        control = .ctl,
        gradient = TRUE,
        scale = "natural"
      ))
      .h0 <- if (.et == "central") (3 * .ef)^(1 / 3) else 2 / sqrt(3) * sqrt(.ef)
      expect_identical(nlmSolveR(.x + 1e-12 * .e), Inf, label = .et)
      .back <- (nlmSolveR(.x) - nlmSolveR(.x - .h0 * .e)) / .h0
      expect_true(is.finite(.back) && .back != 0, label = .et)
      expect_equal(nlmLikEvalC_(.x)$grad[1], .back, tolerance = 1e-6, label = .et)
      .nlmFreeEnv()
    }
  })

  test_that("a theta whose differences solve on neither side fails the evaluation, not a zero slope", {
    skip_on_cran()
    withr::defer(.nlmFreeEnv())
    # the log-likelihood is -Inf at any fdepot but 0.9: the step search forms no
    # difference for fdepot, so nlmixr2NlmEval() reports the point as failed
    # (+Inf, one bad evaluation) instead of a zero fdepot gradient
    .mod <- function() {
      ini({
        fdepot <- 0.9
        tka <- 0.45
        tcl <- 1
        tv <- 3.45
        add.sd <- 0.7
      })
      model({
        ka <- exp(tka)
        cl <- exp(tcl)
        v <- exp(tv)
        d / dt(depot) <- -ka * depot
        f(depot) <- fdepot
        d / dt(center) <- ka * depot - cl / v * center
        cp <- center / v
        ll(err) ~ -0.5 * ((DV - cp) / add.sd)^2 - log(add.sd) - exp(1e300 * (fdepot - 0.9)^2)
      })
    }
    for (.et in c("central", "forward")) {
      .ctl <- nlmControl(print = 0L, eventType = .et, eventSens = "fd")
      .x <- suppressMessages(nlmObjectiveSetup(
        .mod,
        nlmixr2data::theo_sd,
        control = .ctl,
        gradient = TRUE,
        scale = "natural"
      ))
      .e <- c(1, 0, 0, 0, 0)
      expect_true(is.finite(nlmSolveR(.x)), label = .et)
      expect_identical(nlmSolveR(.x + 1e-12 * .e), Inf, label = .et)
      expect_identical(nlmSolveR(.x - 1e-12 * .e), Inf, label = .et)
      .r <- nlmLikEvalC_(.x)
      expect_identical(.r$nBad, 1L, label = .et)
      expect_identical(.r$value, Inf, label = .et)
      .nlmFreeEnv()
    }
  })

  test_that("an eta whose differences solve on neither side keeps the fit and is reported", {
    skip_on_cran()
    # eventSens = "fd" differences eta.f (it enters F); the log-likelihood is -Inf
    # at any eta.f but 0, so no difference forms for it.  Its slope (and its eta
    # Hessian column) are taken as 0, as before, and the fit says how often.
    .mod <- function() {
      ini({
        tka <- 0.45
        tcl <- 1
        tv <- 3.45
        tf <- 0
        eta.cl ~ 0.3
        eta.f ~ 0.1
        add.sd <- 0.7
      })
      model({
        ka <- exp(tka)
        cl <- exp(tcl + eta.cl)
        v <- exp(tv)
        d / dt(depot) <- -ka * depot
        f(depot) <- exp(tf + eta.f)
        d / dt(center) <- ka * depot - cl / v * center
        cp <- center / v
        ll(err) ~ -0.5 * ((DV - cp) / add.sd)^2 - log(add.sd) - exp(1e300 * eta.f^2)
      })
    }
    .fit <- suppressMessages(nlmixr2(
      .mod,
      nlmixr2data::theo_sd,
      "focei",
      foceiControl(print = 0, eventSens = "fd", covMethod = "", maxOuterIterations = 0)
    ))
    expect_true(is.finite(.fit$objective))
    expect_identical(unique(.fit$eta$eta.f), 0)
    .w <- grep("^eta finite differences", .fit$runInfo, value = TRUE)
    expect_length(.w, 1L)
    # 1140 from the eta gradient, 24 from the eta Hessian (12 subjects, twice)
    expect_identical(.w, "eta finite differences solved on neither side 1164 times; slope 0")
  })

  test_that("shi21Central() falls back to a backward difference, and to NA when none forms", {
    .ef <- 7e-7
    .h0 <- (3 * .ef)^(1 / 3)
    # finite only at and below 0: every probe at x + h fails, so the gradient is
    # the backward difference at the starting step, exact for a linear f
    .f <- function(x) if (x[1] > 0) NaN else c(2 * x[1], 3 * x[1])
    .r <- shi21CentralWrap(.f, c(0, 1), .f(c(0, 1)), 1L, .ef)
    expect_equal(drop(.r$gr), c(2, 3), tolerance = 1e-10)
    expect_equal(.r$h, .h0, tolerance = 1e-12)
    # finite at 0 only: no difference forms, and the gradient is NA, not 0
    .f <- function(x) if (x[1] != 0) NaN else c(0, 0)
    .r <- shi21CentralWrap(.f, c(0, 1), .f(c(0, 1)), 1L, .ef)
    expect_identical(drop(.r$gr), c(NA_real_, NA_real_))
    expect_equal(.r$h, .h0, tolerance = 1e-12)
    # finite everywhere: the central difference
    .f <- function(x) c(2 * x[1], 3 * x[1])
    .r <- shi21CentralWrap(.f, c(0, 1), .f(c(0, 1)), 1L, .ef)
    expect_equal(drop(.r$gr), c(2, 3), tolerance = 1e-10)
  })

  # shi21SearchWrap(f, t, f0, idx, ef, forward, maxiter): either step search, with
  # the rl/ru/nu FOCEi's eta search uses
  .ef <- 7e-7
  .h0c <- (3 * .ef)^(1 / 3)
  .h0f <- 2 / sqrt(3) * sqrt(.ef)

  test_that("a backward difference gives way to any other difference that forms", {
    # finite only at and below 0.001: the first steps fail on the plus side, and a
    # shrunk step that solves replaces the backward difference
    .f <- function(x) if (x[1] > 0.001) NaN else 2 * (x[1]^2 + x[1])
    .r <- shi21SearchWrap(.f, 0, .f(0), 1L, .ef, FALSE, 15L)
    expect_equal(.r$h, 0.000249831021, tolerance = 1e-8)
    expect_equal(drop(.r$gr), 2, tolerance = 1e-12) # the central difference is exact
    .r <- shi21SearchWrap(.f, 0, .f(0), 1L, .ef, TRUE, 15L)
    expect_equal(.r$h, 0.000241522946, tolerance = 1e-8)
    expect_equal(drop(.r$gr), 2 * (1 + .r$h), tolerance = 1e-12) # the forward one
  })

  test_that("maxiter = 0 differences at the starting step, without a search", {
    .g <- function(x) 2 * x[1] + 3
    for (.fw in c(FALSE, TRUE)) {
      .r <- shi21SearchWrap(.g, 1, .g(1), 1L, .ef, .fw, 0L)
      expect_equal(.r$h, if (.fw) .h0f else .h0c, tolerance = 1e-12, label = .fw)
      expect_equal(drop(.r$gr), 2, tolerance = 1e-10, label = .fw)
      # the plus side fails: the backward difference
      .g1 <- function(x) if (x[1] > 1) NaN else 2 * x[1]
      expect_equal(drop(shi21SearchWrap(.g1, 1, 2, 1L, .ef, .fw, 0L)$gr), 2, tolerance = 1e-10, label = .fw)
      # neither side: NA
      .g2 <- function(x) if (x[1] != 1) NaN else 0
      expect_identical(drop(shi21SearchWrap(.g2, 1, 0, 1L, .ef, .fw, 0L)$gr), NA_real_, label = .fw)
    }
  })

  test_that("the forward search stops at a step too small to move t on either side", {
    # finite only at t: past about 44 halvings t - h rounds to t at these t
    .only <- function(t0) function(x) if (x[1] != t0) NaN else 0
    for (.t in c(-1, -2, -0.5, -4)) {
      .r <- shi21SearchWrap(.only(.t), .t, 0, 1L, .ef, TRUE, 60L)
      expect_identical(drop(.r$gr), NA_real_, label = .t)
    }
  })

  test_that("a starting step too small to move t may still grow", {
    # at 1e13 the starting steps round to t; the search grows them to hMax = 0.5
    .e <- function(x) exp(x[1] / 1e14)
    .t <- 1e13
    .r <- shi21SearchWrap(.e, .t, .e(.t), 1L, .ef, TRUE, 15L)
    expect_identical(.r$h, 0.5)
    expect_identical(drop(.r$gr), (.e(.t + 0.5) - .e(.t)) / 0.5)
    .r <- shi21SearchWrap(.e, .t, .e(.t), 1L, .ef, FALSE, 15L)
    expect_identical(.r$h, 0.5)
    expect_identical(drop(.r$gr), (.e(.t + 0.5) - .e(.t - 0.5)) / 1)
  })

  # finite only at fdepot = 0.9: fdepot's differences solve on neither side
  .shiEdge2 <- function() {
    ini({
      fdepot <- 0.9
      tka <- 0.45
      tcl <- 1
      tv <- 3.45
      add.sd <- 0.7
    })
    model({
      ka <- exp(tka)
      cl <- exp(tcl)
      v <- exp(tv)
      d / dt(depot) <- -ka * depot
      f(depot) <- fdepot
      d / dt(center) <- ka * depot - cl / v * center
      cp <- center / v
      ll(err) ~ -0.5 * ((DV - cp) / add.sd)^2 - log(add.sd) - exp(1e300 * (fdepot - 0.9)^2)
    })
  }
  test_that("every nlm-family optimizer stops on a gradient whose differences solve on neither side", {
    skip_on_cran()
    .msg <- "the gradient of 'fdepot' is not finite: its finite differences solved on neither side"
    .ctl <- list(
      nlm = nlmControl(print = 0L, eventSens = "fd"),
      nlminb = nlminbControl(print = 0L, eventSens = "fd"),
      optim = optimControl(method = "BFGS", print = 0L, eventSens = "fd"),
      optim = optimControl(method = "CG", print = 0L, eventSens = "fd"),
      n1qn1 = n1qn1Control(print = 0L, eventSens = "fd"),
      lbfgsb3c = lbfgsb3cControl(print = 0L, eventSens = "fd"),
      trust = trustControl(print = 0L, eventSens = "fd")
    )
    for (.i in seq_along(.ctl)) {
      .e <- tryCatch(
        suppressWarnings(suppressMessages(
          nlmixr2(.shiEdge2, nlmixr2data::theo_sd, names(.ctl)[.i], control = .ctl[[.i]])
        )),
        error = conditionMessage
      )
      expect_identical(.e, .msg, label = names(.ctl)[.i])
    }
  })

  test_that("shi21maxFD = 0 differences eta at the starting step, and reports nothing", {
    skip_on_cran()
    .mod <- function() {
      ini({
        tka <- 0.45
        tcl <- 1
        tv <- 3.45
        tf <- 0
        eta.cl ~ 0.3
        eta.f ~ 0.1
        add.sd <- 0.7
      })
      model({
        ka <- exp(tka)
        cl <- exp(tcl + eta.cl)
        v <- exp(tv)
        d / dt(depot) <- -ka * depot
        f(depot) <- exp(tf + eta.f)
        d / dt(center) <- ka * depot - cl / v * center
        cp <- center / v
        cp ~ add(add.sd)
      })
    }
    .fit <- suppressMessages(nlmixr2(
      .mod,
      nlmixr2data::theo_sd,
      "focei",
      foceiControl(print = 0, eventSens = "fd", covMethod = "", maxOuterIterations = 0, shi21maxFD = 0L)
    ))
    expect_true(is.finite(.fit$objective))
    expect_false(any(grepl("^eta finite differences", .fit$runInfo)))
    expect_false(all(.fit$eta$eta.f == 0))
  })

  test_that("innerOpt = 'trust' seeds its quasi-Newton Hessian with a finite matrix", {
    skip_on_cran()
    # as the eta test above; an eta.f Hessian column that solved on neither side
    # was the quasi-Newton seed as NaN, so SR1 skipped every update and Bofill
    # spread it (OFV -565)
    .mod <- function() {
      ini({
        tka <- 0.45
        tcl <- 1
        tv <- 3.45
        tf <- 0
        eta.cl ~ 0.3
        eta.f ~ 0.1
        add.sd <- 0.7
      })
      model({
        ka <- exp(tka)
        cl <- exp(tcl + eta.cl)
        v <- exp(tv)
        d / dt(depot) <- -ka * depot
        f(depot) <- exp(tf + eta.f)
        d / dt(center) <- ka * depot - cl / v * center
        cp <- center / v
        ll(err) ~ -0.5 * ((DV - cp) / add.sd)^2 - log(add.sd) - exp(1e300 * eta.f^2)
      })
    }
    for (.hm in c("bfgs", "sr1", "bofill")) {
      .fit <- suppressMessages(nlmixr2(
        .mod,
        nlmixr2data::theo_sd,
        "focei",
        foceiControl(
          print = 0,
          eventSens = "fd",
          covMethod = "",
          maxOuterIterations = 0,
          innerOpt = "trust",
          hessianMethod = .hm
        )
      ))
      expect_equal(.fit$objective, 122.316911911, tolerance = 1e-9, label = .hm)
      expect_identical(
        grep("^eta finite", .fit$runInfo, value = TRUE),
        "eta finite differences solved on neither side 596 times; slope 0",
        label = .hm
      )
    }
  })
})
