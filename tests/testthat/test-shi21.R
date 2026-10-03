nmTest({
  # shi21Forward()/shi21Central() (src/shi21.cpp): the Shi et al. (2021) step
  # search when one of its outer probes is not finite.

  # eventSens = "fd" finite-differences fdepot (a dosing parameter) for every
  # subject; the log-likelihood is -Inf once fdepot passes `edge` (it starts at
  # 0.9), so a probe past it is not finite.  (rxode2 reads a model from its
  # source when it has one, which still says .(edge).)
  .shiMod <- function(edge) {
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
        ll(err) ~ -0.5 * ((DV - cp) / add.sd)^2 - log(add.sd) - exp(1e9 * (fdepot - .(edge)))
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
    # central: f(x +- h0) are finite and f(x + 3 h0) is not.  The difference
    # was divided by the shrunken step 2 h0 / 3, 1.5 times too large.
    .r <- .shiFirstGrad(0.9 + 0.05, "central")
    expect_false(is.finite(.r[["outer"]]))
    expect_equal(.r[["grad"]], .r[["diff"]], tolerance = 1e-6)
    # forward: f(x + h0) is finite and f(x + 4 h0) is not.  The step grew to
    # 3.5 h0, and the difference at h0 was divided by it (0.28 of the
    # derivative).  Now it is the forward difference at h0, 1.5% off.
    .r <- .shiFirstGrad(0.9 + 0.006, "forward")
    expect_false(is.finite(.r[["outer"]]))
    expect_equal(.r[["grad"]], .r[["fine"]], tolerance = 0.02)
  })

  test_that("the step shi21Central() returns is the one its gradient differences (issue 1140)", {
    # x^2 at 1, finite only on [1, 1.01]: the first probe (1 + h0, h0 = 0.03)
    # fails, the search shrinks to h0 / 6 and takes a forward difference there
    # (the backward probe fails), and every later backward probe fails too.
    # The step returned was h0, at which the function is not finite.
    .f <- function(x) if (x < 1 || x > 1.01) NA_real_ else x^2
    .r <- shi21CentralWrap(.f, 1, 1, 1L, 0.03^3 / 3)
    expect_equal(.r$h, 0.005, tolerance = 1e-12)
    # the forward difference of x^2 at 1 is 2 + h
    expect_equal(drop(.r$gr), 2 + .r$h, tolerance = 1e-12)
  })
})
