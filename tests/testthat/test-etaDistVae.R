## Declared non-normal between-subject distributions under est="vae".
##
## The ELBO does NOT change, and that is the claim under test here.
## `rxEtaDistExpand()` leaves the LATENT standard normal -- `rxz.cl ~ fix(1)`
## has no `theta + eta` form, so `.vaeDataPrep()` marks it `isFree` (zPop 0,
## held there) and reads `omegaFix` from the ini() `fix` -- which is exactly the
## N(0,1) the prior term and the KL in `vaeElboCore` are written for.  The
## non-normality lives in a decoder line inside the inner problem, and rxode2
## differentiates the inverse CDF exactly, so the inner and outer gradients
## already carry `d(eta.declared)/d(eta.latent)`.
##
## Each of those sentences is asserted below rather than argued, because every
## one of them fails SILENTLY: a wrong prior still trains, a missing chain rule
## still returns a finite gradient, and a declared theta that never reaches
## `regressNames` simply keeps its ini() value to the end of the fit.

nmTest({

  ## the declaration + an ODE (linCmt() is out of analytic-gradient scope, which
  ## is a pre-existing limit of .foceiOuterDirs and nothing to do with dist())
  .edVaeMod <- function() {
    function() {
      ini({
        lka <- 0.5; lclm <- 1.5; lv1m <- 1.5
        lclrv <- -1; lv1rv <- -1
        eta.cl + eta.v ~ c(1, 0.3, 1)
        dist(eta.cl) ~ dgamma(shape = 1/exp(lclrv), rate = 1/(exp(lclrv)*exp(lclm)))
        dist(eta.v) ~ dgamma(shape = 1/exp(lv1rv), rate = 1/(exp(lv1rv)*exp(lv1m)))
        eta.ka ~ 0.1
        prop.sd <- 0.3
      })
      model({
        ka <- exp(lka + eta.ka); cl <- eta.cl; v <- eta.v
        d/dt(depot) <- -ka*depot
        d/dt(central) <- ka*depot - (cl/v)*central
        cp <- central/v
        cp ~ prop(prop.sd)
      })
    }
  }

  ## plain Gaussian twin: same structure, mu-referenced.  Used to CALIBRATE the
  ## gradient conventions, so the declared assertions are relative to a known-good
  ## reference rather than to a sign guessed from the source.
  .edVaeNorm <- function() {
    function() {
      ini({
        lka <- 0.5; lcl <- 1.5; lv <- 1.5
        eta.cl + eta.v ~ c(1, 0.3, 1)
        eta.ka ~ 0.1
        prop.sd <- 0.3
      })
      model({
        ka <- exp(lka + eta.ka); cl <- exp(lcl + eta.cl); v <- exp(lv + eta.v)
        d/dt(depot) <- -ka*depot
        d/dt(central) <- ka*depot - (cl/v)*central
        cp <- central/v
        cp ~ prop(prop.sd)
      })
    }
  }

  ## the ui the vae actually sees: expanded, with the stash the expansion destroys
  .edVaeUi <- function(f) {
    .ui <- rxode2::rxUiDecompress(nlmixr2est::nlmixr2(f))
    .st <- .etaDistDeclStash(.ui, rxode2::rxUiEtaDists(.ui))
    .u2 <- rxode2::rxUiDecompress(rxode2::rxEtaDistExpand(.ui))
    .etaDistDeclSet(.u2, .st)
    .u2
  }

  ## ---------------------------------------------------------------- 2.1/2.2 --

  test_that("vae declares etaDist support and emvi/fbvi still do not", {
    expect_true(.etaDistMethodAttr("vae", NULL))
    ## the refusal was never about the ELBO; it is lifted, and only for vae
    expect_false(.etaDistMethodAttr("emvi", NULL))
    expect_false(.etaDistMethodAttr("fbvi", NULL))
    ## the nonparametric methods stay refused -- a declared family contradicts them
    expect_false(.etaDistMethodAttr("npag", NULL))
  })

  ## ------------------------------------------------------------------- 2.3 --

  test_that("vaeControl carries etaDistWarmStart", {
    ## Real gap, not bookkeeping: `.preProcessEtaDist` treats ABSENCE as TRUE, so
    ## before this every declared vae fit would have run a nested est="saem"
    ## surrogate with no way to turn it off.
    expect_true(vaeControl()$etaDistWarmStart)
    expect_false(vaeControl(etaDistWarmStart = FALSE)$etaDistWarmStart)
    expect_error(vaeControl(etaDistWarmStart = "yes"))
    expect_error(vaeControl(etaDistWarmStart = c(TRUE, TRUE)))
    ## and it round-trips.  UNFILTERED on purpose: getValidNlmixrCtl.vae does
    ## `do.call(vaeControl, .ctl)` on the whole returned list (R/vae.R:747), so
    ## every element must be a declared argument.  Filtering to
    ## names(formals(vaeControl)) first would drop exactly the elements that
    ## would throw, leaving a test that cannot fail for the reason it exists.
    .c <- vaeControl(etaDistWarmStart = FALSE)
    expect_false(do.call(vaeControl, .c)$etaDistWarmStart)
  })

  ## ---------------------------------------------------------------- Risk 1 --

  test_that(".preProcessEtaDist runs first in the hook chain", {
    ## Everything below depends on it.  If it does not run first, `.vaeNonMuThetas`
    ## sees the UNEXPANDED model, the declared thetas never reach `regressNames`,
    ## and the fit completes with the declaration frozen at its ini() values -- no
    ## error, no warning, wrong answer.
    .h <- .orderPreProcessHookNames(c(".preProcessBoundedTransform", ".preProcessIov",
                                      ".preProcessEtaDist", ".preProcessVaeNonMuTheta"))
    expect_identical(.h[1], ".preProcessEtaDist")
    expect_identical(.h[length(.h)], ".preProcessBoundedTransform")
    ## and in the registered chain itself
    .r <- .orderPreProcessHookNames(ls(.preProcessHooks))
    if (".preProcessEtaDist" %in% .r) expect_identical(.r[1], ".preProcessEtaDist")
  })

  ## -------------------------------------------------- the ELBO prior terms --

  test_that("the expansion leaves the vae latent standard normal", {
    skip_on_cran()
    .u <- .edVaeUi(.edVaeMod())
    .p <- .vaeDataPrep(.u, nlmixr2data::theo_sd, vaeControl(covariateSelection = FALSE))
    .z <- grep("^rxz[.]", .p$etaNames)
    expect_gt(length(.z), 0L)                      # the declaration produced latents
    ## isFree: no `theta + eta` form, so zPop is forced to 0 and HELD there
    expect_true(all(.p$isFree[.z]))
    expect_equal(unname(.p$zPop[.z]), rep(0, length(.z)))
    ## unit variance, and FIXED -- this is the N(0,1) the KL is written for
    expect_equal(unname(.p$omega[.z]), rep(1, length(.z)))
    expect_true(all(.p$omegaFix[.z]))
    ## the undeclared eta is untouched by any of this
    .k <- match("eta.ka", .p$etaNames)
    expect_false(is.na(.k))
    expect_false(.p$omegaFix[.k])
  })

  test_that("training leaves the declared latent's prior alone", {
    skip_on_cran()
    ## A fixed omega that the M-step quietly estimates anyway would make the KL
    ## and the prior disagree, and nothing else in the fit would say so.
    r <- suppressWarnings(suppressMessages(
      nlmixr2(.edVaeMod()(), nlmixr2data::theo_sd, est = "vae",
              control = vaeControl(print = 0L, calcTables = FALSE, returnVae = TRUE,
                                   covariateSelection = FALSE, etaDistWarmStart = FALSE,
                                   itersBurnIn = 10L, iters = 20L, klWarmup = 5L,
                                   gammaIter = 15L))))
    ## returnVae=TRUE short-circuits .vaeToFit() and hands back .vaeTrain()'s raw
    ## list, where `omega` is an unnamed numeric vector -- `omegaMat` is the
    ## dimnamed matrix (R/vaeFit.R:320).  rownames(a vector) is NULL, so reading
    ## `omega` here would have made every assertion below vacuous.
    .om <- r$omegaMat
    .z <- grep("^rxz[.]", rownames(.om), value = TRUE)
    expect_gt(length(.z), 0L)
    for (.n in .z) expect_equal(unname(.om[.n, .n]), 1, tolerance = 1e-8)
    ## The expansion DROPS the block off-diagonals -- the correlation moved into
    ## the `rxCor.*` thetas.  If the omega M-step filled one back in, the copula
    ## would have two owners and the fit would be quietly over-parameterized.
    if (length(.z) > 1L) {
      for (.i in seq_along(.z)) for (.j in seq_along(.z)) if (.i != .j) {
        expect_equal(unname(.om[.z[.i], .z[.j]]), 0, tolerance = 1e-8)
      }
    }
  })

  ## ------------------------------------------------------ the regressed set --

  test_that("the declared thetas and rxCor.* reach the regressed set", {
    skip_on_cran()
    ## If they do not, the M-step never touches them and the fit returns ini().
    .u <- .edVaeUi(.edVaeMod())
    .nm <- .vaeNonMuThetas(.u)
    expect_true(all(c("lclm", "lclrv", "lv1m", "lv1rv") %in% .nm))
    ## rxCor.* is the Gaussian copula parameter; vae regresses it, and if it comes
    ## back NA the whole M-step declines (R/vaeGrad.R)
    expect_gt(length(grep("^rxCor[.]", .nm)), 0L)
    ## lka stays OUT -- it is mu-referenced through eta.ka
    expect_false("lka" %in% .nm)
    ## and .vaeDataPrep agrees with the helper
    .p <- .vaeDataPrep(.u, nlmixr2data::theo_sd, vaeControl(covariateSelection = FALSE))
    expect_true(all(.nm %in% .p$regressNames))
  })

  ## ------------------------------------------------------------------- 2.4 --

  test_that("the cheap scope probe accepts a declared model", {
    skip_on_cran()
    ## `.vaeGradInScope()` is a direction-set probe with no symengine/gcc pass, so
    ## what it answers is "is this model SHAPED like one the analytic gradient
    ## handles" -- and a declared model is.  Whether the augmented model actually
    ## BUILDS is a separate question, answered later by `ui$foceiOuter`, and today
    ## it does NOT build for a declared model.  So the safety property is not
    ## enforced here; it is enforced on the real build result -- see the segfault
    ## regression test below.
    expect_true(.vaeGradInScope(.edVaeUi(.edVaeMod())))
    expect_true(.vaeGradInScope(rxode2::assertRxUi(.edVaeNorm()())))
  })

  test_that("grad without a built augmented model downgrades instead of crashing", {
    skip_on_cran()
    ## REGRESSION TEST for a segfault.  With nonMuTheta="grad" and no augmented
    ## model registered, nothing sized the shared solve pool for it, and the first
    ## burn-in ELBO step died in iniSubject() -> _setIndPointersByThread() on a
    ## null gInfusionRate -- one solve after vaeInnerUpdatePar_.
    ##
    ## A declared model is the case that REACHES this (`ui$foceiOuter` does not
    ## build for one), but the bug was never about declarations: ANY model whose
    ## augmented build declines under "grad" hit it.  The fix downgrades on the
    ## actual build result in `.vaeFitModel`, which owns the only control that
    ## reaches both `.vaeGradInit` and `vaeTrainCpp_`'s `useGrad` gate.
    r <- suppressWarnings(suppressMessages(
      nlmixr2(.edVaeMod()(), nlmixr2data::theo_sd, est = "vae",
              control = vaeControl(nonMuTheta = "grad", print = 0L,
                                   calcTables = FALSE, returnVae = TRUE,
                                   covariateSelection = FALSE,
                                   etaDistWarmStart = FALSE, itersBurnIn = 10L,
                                   iters = 20L, klWarmup = 5L, gammaIter = 15L))))
    ## it COMPLETED -- that is the property under test
    expect_false(is.null(r))
    ## and it completed on the regression, rather than pretending grad ran
    expect_equal(r$nRegGrad, 0L)
    expect_equal(r$nRegFallback, 0L)
    expect_true(all(is.finite(r$regressTheta)))
  })

  ## ------------------------------------------------- inner gradient vs c.d. --

  test_that("the inner eta-gradient carries the inverse-CDF chain rule", {
    skip_on_cran()
    ## THE substantive claim of declared-distribution support on the FOCEi-family
    ## inner problem: d(eta.declared)/d(eta.latent) through
    ## gammapInv(shape, phiU(latent)) must reach lpInner.  rxode2 differentiates
    ## the inverse CDF exactly (.rxD$gammapInv gives dq/dp = 1/gammapDer(...)),
    ## so this should hold to solver accuracy, not merely approximately.
    ##
    ## foceiInnerLp is d(likInner)/d(eta) up to one fixed convention.  Calibrate
    ## it on the Gaussian twin, then require the SAME convention on the declared
    ## model: what is under test is the chain rule, not the sign of lpInner.
    ##
    ## MEASURED tolerances, and they are not the ones you would guess -- the
    ## DECLARED columns are the accurate ones (ratio 1.0000, agreeing to ~1e-6),
    ## while the ordinary mu-referenced eta.ka column is the loose one at ~1.4e-3
    ## because exp(lka + eta.ka) feeds d/dt and carries the ODE solver's own
    ## error.  Asserting one pooled tolerance would either fail on eta.ka or go
    ## slack enough to miss a broken chain rule.
    .cd <- function(id, eta, h) {
      vapply(seq_along(eta), function(j) {
        .hj <- h * max(1, abs(eta[j]))
        .up <- eta; .up[j] <- .up[j] + .hj
        .dn <- eta; .dn[j] <- .dn[j] - .hj
        (likInner(.up, id) - likInner(.dn, id)) / (2 * .hj)
      }, numeric(1))
    }
    .ctl <- vaeControl(covariateSelection = FALSE)
    .d <- nlmixr2data::theo_sd
    .N <- length(unique(.d$ID))

    ## --- calibration on the Gaussian twin -------------------------------------
    .un <- rxode2::assertRxUi(.edVaeNorm()())
    .e0 <- matrix(0.05, .N, 3L)
    .ienv <- .vaeInnerSetup(.un, .d, .e0, .ctl)
    .lpN <- as.numeric(foceiInnerLp(.e0[1, ], 1L))
    .cdN <- .cd(1L, .e0[1, ], 1e-4)
    .vaeInnerFree()
    .k <- stats::median(.lpN / .cdN)
    ## lpInner IS d(likInner)/d(eta): same sign, unit scale, no prior offset
    expect_equal(.k, 1, tolerance = 1e-3)

    ## --- the same convention on the declared model ----------------------------
    ## gammapInv is an ITERATIVE inversion, so its own solve tolerance floors any
    ## finite difference through it; a relative step of 1e-4 clears that and
    ## stays in the linear regime (measured identical at 1e-5 and 1e-6).
    .ud <- .edVaeUi(.edVaeMod())
    .p <- .vaeDataPrep(.ud, .d, .ctl)
    .e1 <- matrix(0.05, .N, .p$zDim)
    .ienv <- .vaeInnerSetup(.ud, .d, .e1, .ctl)
    on.exit(.vaeInnerFree(), add = TRUE)
    .lpD <- as.numeric(foceiInnerLp(.e1[1, ], 1L))
    .cdD <- .cd(1L, .e1[1, ], 1e-4)
    expect_true(all(is.finite(.lpD)))
    .zi <- grep("^rxz[.]", .p$etaNames)
    .oi <- setdiff(seq_along(.p$etaNames), .zi)
    expect_gt(length(.zi), 0L); expect_gt(length(.oi), 0L)
    ## Declared and undeclared asserted SEPARATELY: pooled, a dead chain rule on
    ## the two declared latents can hide behind a correct eta.ka.
    expect_equal(.lpD[.zi], .k * .cdD[.zi], tolerance = 1e-3)
    expect_equal(.lpD[.oi], .k * .cdD[.oi], tolerance = 1e-2)
    ## and the declared columns must not be ZERO -- a missing dq/dp reads as "no
    ## dependence", which would pass any tolerance against its own zero reference
    expect_gt(min(abs(.lpD[.zi])), 1e-8)
  })

  ## ------------------------------------------------- outer gradient vs c.d. --

  ## No outer-gradient test here.  `.vaeGradInScope()` now declines a declared
  ## model (it segfaults -- see the note there), so `.vaeGradInit`/`.vaeGradEval`
  ## on one is a route no fit can take, and asserting against it would pin
  ## behavior nothing depends on.  When the crash is fixed, the test to add is
  ## that `.g` is non-NULL and covers every name in `regressNames` INCLUDING
  ## `rxCor.*`: `.vaeGradEval` subsets by name and returns NULL on any NA, so a
  ## copula theta the direction set does not carry silently disables the whole
  ## M-step.

  ## ------------------------------------------------------------ end to end --

  test_that("a declared vae fit runs end to end and its thetas move", {
    skip_on_cran()
    ## Plausible estimates are not the assertion: a silent downgrade to bobyqa,
    ## or a declaration frozen at ini(), both still produce a finished fit.
    ## Default nonMuTheta ("regress"), which is the route a declared model
    ## actually takes -- "grad" is refused for one (see the downgrade test).
    r <- suppressWarnings(suppressMessages(
      nlmixr2(.edVaeMod()(), nlmixr2data::theo_sd, est = "vae",
              control = vaeControl(print = 0L, calcTables = FALSE,
                                   returnVae = TRUE, covariateSelection = FALSE,
                                   etaDistWarmStart = FALSE,
                                   itersBurnIn = 10L, iters = 30L, klWarmup = 5L,
                                   gammaIter = 20L))))
    ## The bobyqa regression owns these thetas on this route, so the analytic
    ## gradient must NOT have fired -- and it must not have "fallen back" either,
    ## which would mean it was attempted.
    expect_equal(r$nRegGrad, 0L)
    expect_equal(r$nRegFallback, 0L)
    ## And the declared thetas MOVED off their ini() values.  Read from
    ## `regressTheta`, which is setNames()d over prep$regressNames
    ## (R/vaeFit.R:328) -- `r$theta` is NULL on this object, and NULL[["lclm"]]
    ## returns NULL rather than erroring, so asserting through it would have
    ## PASSED on a declaration frozen at ini(): exactly the failure this file's
    ## header warns about.
    .rt <- r$regressTheta
    expect_false(is.null(.rt))
    expect_true(all(c("lclm", "lclrv", "lv1m", "lv1rv") %in% names(.rt)))
    expect_true(all(is.finite(.rt[c("lclm", "lclrv", "lv1m", "lv1rv")])))
    expect_false(isTRUE(all.equal(unname(.rt[["lclm"]]), 1.5)))
  })

})
