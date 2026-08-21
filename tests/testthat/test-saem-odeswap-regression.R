nmTest({
  # Phase 4 (see the SAEM general-likelihood theta plan) converts SAEM's own
  # per-individual solve/read call sites (`setupRx`/`user_function`,
  # src/saem.cpp) to route through the shared `odeSwap` pool -- but ONLY for a
  # general-likelihood fit (distribution==4); every normal/poisson/binomial
  # SAEM fit's own setup/solve code path is meant to stay completely
  # untouched. This test pins a normal-error SAEM fit's exact numeric result
  # (single-threaded, fixed seed, so it is genuinely deterministic -- not just
  # "close") BEFORE that conversion, as the "zero blast radius for existing
  # fits" proof the plan's evaluation criteria require (item 10): if a later
  # change to setupRx/user_function's shared code paths (not gated on
  # distribution==4) ever perturbs a normal fit's result, this goes red.
  test_that("normal-error SAEM fit is unaffected by the general-lik odeSwap plumbing", {
    rxode2::setRxThreads(1)
    mod <- function() {
      ini({
        tka <- 0.45; tcl <- 1; tv <- 3.45
        eta.ka ~ 0.6; eta.cl ~ 0.3; eta.v ~ 0.1
        add.sd <- 0.7
      })
      model({
        ka <- exp(tka + eta.ka); cl <- exp(tcl + eta.cl); v <- exp(tv + eta.v)
        d / dt(depot) <- -ka * depot
        d / dt(center) <- ka * depot - cl / v * center
        cp <- center / v
        cp ~ add(add.sd)
      })
    }
    .ctl <- saemControl(nBurn = 10, nEm = 10, print = 0, seed = 42)
    f <- suppressWarnings(suppressMessages(
      .nlmixr(mod, nlmixr2data::theo_sd, est = "saem", control = .ctl)
    ))
    expect_equal(f$objf, 114.912042038049947, tolerance = 1e-8)
    .fx <- unname(fixef(f))
    expect_equal(.fx,
                 c(0.470404856283517, 1.005906363801880,
                   3.456526558949254, 0.701092895255751),
                 tolerance = 1e-8)
  })
})
