## Slow (weekly batch): full grad-vs-regress fits on theo_sd.
##
## Both modes drive the SAME (full outer) objective -- the Laplace determinant +
## 0.5*log|Omega^-1| + the DV-transform Jacobian -- with every mu-referenced theta
## held at its current M-step value.  "regress" fully re-solves the regressed
## thetas with a bounded bobyqa search each M-step; "grad" takes one Adam step
## along the exact analytic gradient.  The invariant is that grad reaches an
## objective at least as good as regress: they optimize the same thing, and grad
## has the exact derivative.  Distance to the FOCEi MLE is a sanity check, not
## the invariant -- the VAE objective is not the FOCEi objective, so either mode
## can sit marginally closer to the FOCEi optimum.
##
## The residual (add.sd) is asserted on purpose.  An error parameter's live value
## is the `a` vector, NOT its theta slot (vaeBuildTh rebuilds that slot from `a`
## on every evaluation), so a grad update written only to the theta slot is
## silently discarded: the fit still converges and reports a plausible tv while
## add.sd and the objective are badly wrong.  A tv-only assertion cannot see it.

nmTest({
  .mod <- function() {
    ini({ tka <- 0.45; tcl <- 1; tv <- c(2, 3.45, 5); add.sd <- c(0, 0.7, 5)
      eta.ka ~ 0.6; eta.cl ~ 0.3 })
    model({ ka <- exp(tka + eta.ka); cl <- exp(tcl + eta.cl); v <- exp(tv)
      d / dt(depot) <- -ka * depot
      d / dt(center) <- ka * depot - cl / v * center
      cp <- center / v
      cp ~ add(add.sd) })
  }
  .ctl <- function(m, ...) {
    vaeControl(nonMuTheta = m, print = 0L, calcTables = FALSE, ...)
  }

  test_that("nonMuTheta='grad' runs the analytic-gradient M-step and fits near the FOCEi MLE", {
    skip_on_cran()
    .mle <- suppressMessages(
      nlmixr2(.mod(), nlmixr2data::theo_sd, est = "focei",
              control = foceiControl(print = 0L, covMethod = "", calcTables = FALSE)))$theta[["tv"]]

    .reg <- suppressWarnings(suppressMessages(rxode2::rxWithSeed(42,
      nlmixr2(.mod(), nlmixr2data::theo_sd, est = "vae",
              control = .ctl("regress", returnVae = TRUE)))))
    .grd <- suppressWarnings(suppressMessages(rxode2::rxWithSeed(42,
      nlmixr2(.mod(), nlmixr2data::theo_sd, est = "vae",
              control = .ctl("grad", returnVae = TRUE)))))

    ## the mechanism actually ran (a silent bobyqa fallback would still produce a
    ## plausible number, so assert the path, not just the value)
    expect_equal(.reg$nRegGrad, 0L)
    expect_true(.grd$nRegGrad > 0L)
    expect_equal(.grd$nRegFallback, 0L)

    ## grad lands close to the FOCEi MLE (measured ~4e-4 here)
    expect_true(abs(.grd$regressTheta[["tv"]] - .mle) < 0.001)
    expect_true(abs(.reg$regressTheta[["tv"]] - .mle) < 0.005)
  })

  test_that("nonMuTheta='grad' estimates the residual and beats 'regress' on the objective", {
    skip_on_cran()
    ## Regression guard for the discarded-error-parameter bug: when the grad
    ## M-step failed to write an error parameter back to `a`, add.sd converged to
    ## ~1.70 against ~0.80 for "regress" and the objective was ~86 units worse,
    ## while tv stayed within 3e-3 of the MLE -- invisible to a tv-only check.
    .fitFor <- function(m) {
      suppressWarnings(suppressMessages(rxode2::rxWithSeed(42,
        nlmixr2(.mod(), nlmixr2data::theo_sd, est = "vae", control = .ctl(m)))))
    }
    .reg <- .fitFor("regress")
    .grd <- .fitFor("grad")

    ## the residual is actually estimated: both modes agree (measured ~0.004 apart)
    expect_equal(.grd$theta[["add.sd"]], .reg$theta[["add.sd"]], tolerance = 0.05)
    ## same objective, exact gradient: grad must not do worse than the
    ## derivative-free search (measured: grad ~0.22 better)
    expect_lte(.grd$objf, .reg$objf + 1)
  })
})
