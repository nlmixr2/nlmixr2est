## `est="vae"` REFUSES a declared non-normal between-subject distribution.
##
## This file used to assert the opposite, on a structural argument that is still
## true sentence by sentence: `rxEtaDistExpand()` leaves the LATENT standard
## normal, which is exactly what `vaeElboCore`'s prior term and KL are written
## for, and rxode2 differentiates the inverse CDF exactly, so the inner and
## outer gradients carry `d(eta.declared)/d(eta.latent)`.  Every one of those
## held.  The fits were still wrong.
##
## Measured on Bauer's four arms (inst/sim/benchArms-results.md):
##
##   arm   CV     MARE mean   MARE rv   rho (truth 0.5)
##   g1    30%      127.8%      99.8%    0.998
##   g3    71%       37.7%     100.0%    0.998
##   g2   100%       52.7%      90.2%   -0.999
##   g4   141%      100.0%     662.6%   -0.999
##
## vae is the ONLY method whose failure does not improve as the declared
## dispersion rises.  The copula is pinned at |rho| = 1 on every arm and
## sign-flips between g3 and g2; the declared relative variance goes to ~1e-4
## (rvCL 9.96e-05 on g1) while the residual inflates to 3-4x truth to absorb it.
## The ELBO relocates the between-subject variability rather than estimating it.
##
## The lesson the removal is worth keeping for: a structural argument that every
## component is correct is not evidence that the fit is.  Each sentence above
## was asserted by a test in this file, every one passed, and the arms failed
## anyway.  Re-enable only behind a fit that recovers them.

nmTest({

  test_that("vae refuses a declared distribution, and emvi/fbvi/npag still do", {
    expect_false(.etaDistMethodAttr("vae", NULL))
    expect_false(.etaDistMethodAttr("emvi", NULL))
    expect_false(.etaDistMethodAttr("fbvi", NULL))
    ## the nonparametric methods stay refused for their own reason -- a declared
    ## family contradicts modelling the distribution nonparametrically
    expect_false(.etaDistMethodAttr("npag", NULL))
    ## and the methods that DO carry it are unaffected by this removal
    expect_true(.etaDistMethodAttr("saem", NULL))
    expect_true(.etaDistMethodAttr("focei", NULL))
  })

  test_that("a declared model under est=\"vae\" is refused by name", {
    .d <- data.frame(name = "eta.cl", etaDist = "dgamma(1, 1)",
                     stringsAsFactors = FALSE)
    expect_error(.etaDistRefuse(.d, "vae", NULL), "vae")
    expect_error(.etaDistRefuse(.d, "vae", NULL), "eta.cl")
  })

  test_that("vaeControl still carries etaDistWarmStart", {
    ## Kept although vae no longer accepts a declaration: it is a declared
    ## argument of vaeControl(), and getValidNlmixrCtl.vae round-trips the WHOLE
    ## returned list through do.call(vaeControl, .ctl) (R/vae.R), so removing it
    ## would break that round trip rather than tidy anything.
    expect_true(vaeControl()$etaDistWarmStart)
    expect_false(vaeControl(etaDistWarmStart = FALSE)$etaDistWarmStart)
    expect_error(vaeControl(etaDistWarmStart = "yes"))
    .c <- vaeControl(etaDistWarmStart = FALSE)
    expect_false(do.call(vaeControl, .c)$etaDistWarmStart)
  })

  test_that(".preProcessEtaDist runs first in the hook chain", {
    ## NOT vae-specific and kept for that reason: every method that carries a
    ## declaration depends on it.  If it does not run first the declared thetas
    ## never reach the regressed set and the fit completes with them frozen at
    ## their ini() values -- no error, no warning, wrong answer.
    .h <- .orderPreProcessHookNames(c(".preProcessBoundedTransform", ".preProcessIov",
                                      ".preProcessEtaDist", ".preProcessVaeNonMuTheta"))
    expect_identical(.h[1], ".preProcessEtaDist")
    expect_identical(.h[length(.h)], ".preProcessBoundedTransform")
    .r <- .orderPreProcessHookNames(ls(.preProcessHooks))
    if (".preProcessEtaDist" %in% .r) expect_identical(.r[1], ".preProcessEtaDist")
  })
})
