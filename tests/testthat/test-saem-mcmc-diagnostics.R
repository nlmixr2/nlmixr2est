# saem computed its MCMC acceptance rate, adapted nothing with it, and threw it
# away.  A chain that had stopped moving therefore looked exactly like one
# exploring properly, and "the chain is not mixing" is a diagnosis that has to
# be made from the fit rather than guessed at.
#
# These are diagnostics only: they must not change a single fitted number.
nmTest({
  .diagMod <- function() {
    ini({
      tka <- 0.45; tcl <- 1; tv <- 3.45
      eta.ka ~ 0.6; eta.cl ~ 0.3; eta.v ~ 0.1
      add.sd <- 0.7
    })
    model({
      ka <- exp(tka + eta.ka); cl <- exp(tcl + eta.cl); v <- exp(tv + eta.v)
      linCmt() ~ add(add.sd)
    })
  }

  test_that("the mixing traces are recorded, named and shaped", {
    skip_on_cran()
    .f <- suppressMessages(nlmixr2(
      .diagMod(),
      theo_sd,
      est = "saem",
      control = saemControl(print = 0, nBurn = 20, nEm = 20, seed = 42L, covMethod = "")
    ))
    .n <- 40L
    # one row per iteration, one column per kernel
    expect_equal(dim(.f$mcmcAccept), c(.n, 3L))
    expect_equal(colnames(.f$mcmcAccept), c("prior", "rw", "coord"))
    expect_true(all(is.finite(.f$mcmcAccept)))
    expect_true(all(.f$mcmcAccept >= 0 & .f$mcmcAccept <= 1))

    # the per-parameter traces carry the sampled parameters' own names, which
    # is the whole point of them: a positional column is not actionable
    # `$` dispatches to the fit environment; `[[` indexes the underlying data
    # frame and would quietly return NULL for every one of these
    for (.m in list(.f$mcmcPhiSd, .f$mcmcPhiAcf, .f$mcmcAcceptCol)) {
      expect_equal(nrow(.m), .n)
      expect_equal(colnames(.m), c("tka", "tcl", "tv"))
    }
    # a lag-1 autocorrelation is a correlation
    expect_true(all(abs(.f$mcmcPhiAcf[-1, ]) <= 1))
    # ...and the stuck fraction is a fraction
    expect_equal(length(.f$mcmcStuck), .n)
    expect_true(all(.f$mcmcStuck >= 0 & .f$mcmcStuck <= 1))
  })

  test_that("kernel 3's acceptance is reported per column, not pooled", {
    # Kernel 3 is Metropolis-within-Gibbs: it proposes ONE coordinate at a
    # time, so its acceptance is already a per-column quantity and the pooled
    # rate averages it away.  That matters -- a fit whose pooled `coord` rate
    # sits exactly on its target can still have one column accepting nearly
    # everything, which is a coordinate whose proposals the likelihood is not
    # rejecting, and whose chain is then exploring the prior.
    skip_on_cran()
    .f <- suppressMessages(nlmixr2(
      .diagMod(),
      theo_sd,
      est = "saem",
      control = saemControl(print = 0, nBurn = 20, nEm = 20, seed = 42L, covMethod = "")
    ))
    .col <- .f$mcmcAcceptCol
    expect_true(all(is.finite(.col)))
    expect_true(all(.col >= 0 & .col <= 1))
    # the pooled `coord` rate is a weighted mean of the columns, so it must lie
    # within their range -- this is what makes the pooled number misleading
    .last <- nrow(.col)
    expect_gte(.f$mcmcAccept[.last, "coord"], min(.col[.last, ]) - 1e-8)
    expect_lte(.f$mcmcAccept[.last, "coord"], max(.col[.last, ]) + 1e-8)
  })

  test_that("recording the diagnostics does not change the fit", {
    # The gate.  These accumulate counters and read phiM; they must not touch a
    # proposal, an acceptance test or an M-step.  Pinned rather than compared
    # to a second run, so a change that perturbed the sampler would show up
    # here even if it were self-consistent.
    skip_on_cran()
    .old <- rxode2::getRxThreads(verbose = FALSE)
    on.exit(rxode2::setRxThreads(.old))
    rxode2::setRxThreads(1L)
    .m <- function() {
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
    .f <- suppressWarnings(suppressMessages(
      nlmixr2(
        .m,
        nlmixr2data::theo_sd,
        est = "saem",
        control = saemControl(print = 0, nBurn = 10, nEm = 10, seed = 42L, calcTables = FALSE, covMethod = "")
      )
    ))
    expect_equal(.f$objf, 115.036204209672, tolerance = 1e-8)
    expect_equal(
      unname(fixef(.f)),
      c(0.454331063367426, 1.011629400082714, 3.456416552020142, 0.702842211703013),
      tolerance = 1e-8
    )
  })
})
