## The DIRECT parameterization: the declared eta is itself the random effect and
## carries its family as a prior, rather than being a standard normal that a
## decoder turns into one.
##
## Every assertion here is written so it FAILS if the Gaussian prior is still in
## use.  That matters more than usual: the broken version does not error, does
## not warn, and returns estimates that look entirely reasonable -- the omega
## entry is a fixed placeholder 1, so a fit reads it as a unit variance and
## completes.  "It ran and the numbers look plausible" is precisely the
## signature of the bug, so no test here is allowed to rest on it.

.edDirectModel <- function() {
  function() {
    ini({
      lclm <- log(4)
      lv <- log(50)
      lclrv <- log(0.8)
      prop.sd <- c(0, 0.2)
      dist(eta.cl) ~ dgamma(shape = 1/exp(lclrv),
                            rate = 1/(exp(lclrv)*exp(lclm)))
    })
    model({
      cl <- eta.cl
      v <- exp(lv)
      d/dt(centr) <- -cl/v*centr
      cp <- centr/v
      cp ~ prop(prop.sd)
    })
  }
}

.edDirectData <- function(n = 40, seed = 42) {
  set.seed(seed)
  .shape <- 2
  .obs <- do.call(rbind, lapply(seq_len(n), function(.i) {
    .cl <- stats::rgamma(1, .shape, .shape/5.104)
    .t <- c(0.25, 0.5, 1, 2, 4, 6, 8, 12, 24)
    data.frame(ID = .i, TIME = .t,
               DV = 100/50*exp(-.cl/50*.t)*exp(stats::rnorm(length(.t), 0, 0.15)),
               AMT = 0, EVID = 0)
  }))
  .d <- rbind(data.frame(ID = seq_len(n), TIME = 0, DV = 0, AMT = 100, EVID = 1),
              .obs)
  .d[order(.d$ID, .d$TIME, -.d$EVID), ]
}

test_that("the route survives to the estimator, on the stash and not on the ui", {
  .f <- .edDirectModel()
  .u0 <- rxode2::as.rxUi(.f)
  for (.p in c("cdf", "direct")) {
    .ui <- .preProcessEtaDist(.u0, est = "saem",
                              control = saemControl(etaDistParam = .p,
                                                    etaDistMstep = TRUE,
                                     etaDistWarmStart = FALSE))$ui
    .ui <- rxode2::rxUiDecompress(.ui)
    expect_equal(.etaDistIsDirect(.ui), identical(.p, "direct"))
    expect_equal(.etaDistEtaPrefix(.ui), if (.p == "direct") "rxd." else "rxz.")
    ## the eta really was renamed the way the prefix says
    .en <- .ui$iniDf$name[!is.na(.ui$iniDf$neta1) &
                            .ui$iniDf$neta1 == .ui$iniDf$neta2]
    expect_true(all(grepl(paste0("^", .etaDistEtaPrefix(.ui)), .en, fixed = FALSE)))
    ## ...and it is the STASH that says so.  `etaDistInfo` records it too, but
    ## that is a ui environment variable and it does not reach the estimator:
    ## before this was moved onto the stash, saem.R read the route as "cdf" off
    ## a ui whose eta was already named `rxd.eta.cl`, looked up `rxz.eta.cl`,
    ## found nothing, and silently fitted a standard normal.
    .st <- .etaDistDeclGet(.ui)
    expect_false(is.null(.st))
    expect_equal(.st$param, .p)
  }
})

test_that("the declared distribution resolves to saem's parameters on both routes", {
  .f <- .edDirectModel()
  .u0 <- rxode2::as.rxUi(.f)
  for (.p in c("cdf", "direct")) {
    .ui <- rxode2::rxUiDecompress(
      .preProcessEtaDist(.u0, est = "saem",
                         control = saemControl(etaDistParam = .p,
                                               etaDistMstep = TRUE,
                                     etaDistWarmStart = FALSE))$ui)
    .en <- .ui$iniDf$name[!is.na(.ui$iniDf$neta1) &
                            .ui$iniDf$neta1 == .ui$iniDf$neta2]
    .edi <- .etaDistMstepInfo(.ui, .ui$saemEtaTrans, .en, .ui$saemParamsToEstimate)
    ## a NULL here is what "no declarations" looks like to the C++ side, and on
    ## the direct route that means the prior itself went missing
    expect_false(is.null(.edi))
    expect_equal(.edi$direct, if (.p == "direct") 1L else 0L)
    expect_equal(length(.edi$latent), 1L)
    expect_true(.edi$latent >= 0)
    expect_equal(.edi$fam, .etaDistFamilyCode("dgamma(2, 2)"))
  }
})

test_that("a declared gamma sampled on the direct route cannot go negative", {
  ## THE engagement test.  A gamma random effect has support (0, Inf), so every
  ## fitted eta must be positive -- kernel 1 draws from the family itself and
  ## kernels 2 and 3 walk on log(eta), so nothing can propose outside it.  Under
  ## the Gaussian prior the same fit produces negative etas, because the omega
  ## placeholder makes it a standard normal.
  ##
  ## Measured on this model at the time it was written:
  ##   Gaussian prior leaking in   eta in [-0.199,  1.977]
  ##   direct prior engaged        eta in [ 0.603, 13.825]
  ##
  ## The bound, not the estimates: an estimate can look reasonable either way.
  skip_on_cran()
  .fit <- suppressWarnings(nlmixr2(
    .edDirectModel(), .edDirectData(), "saem",
    saemControl(nBurn = 20, nEm = 20, nmc = 3, print = 0, seed = 99,
                etaDistParam = "direct", etaDistMstep = TRUE)))
  .eta <- .fit$eta[[2]]
  expect_true(all(.eta > 0),
              info = paste0("declared gamma eta went to ", min(.eta),
                            "; the Gaussian prior is still in use"))
  ## and it is not merely positive by accident of a tiny spread
  expect_gt(max(.eta), 2)
})

test_that("focei and imp refuse the direct route, and say why", {
  .f <- .edDirectModel()
  .u0 <- rxode2::as.rxUi(.f)
  for (.e in c("focei", "imp")) {
    .ctl <- if (.e == "focei") foceiControl(etaDistParam = "direct", etaDistWarmStart = FALSE)
            else impmapControl(etaDistParam = "direct", etaDistWarmStart = FALSE)
    expect_error(.preProcessEtaDist(.u0, est = .e, control = .ctl),
                 "no interior mode")
  }
  ## and saem does NOT refuse it
  expect_error(.preProcessEtaDist(.u0, est = "saem",
                                  control = saemControl(etaDistParam = "direct",
                                                        etaDistMstep = TRUE,
                                     etaDistWarmStart = FALSE)),
               NA)
})

test_that("etaDistMstep does not apply to the direct route", {
  ## This test used to assert a REFUSAL ("direct requires etaDistMstep=TRUE").
  ## That refusal named the wrong cause and has been removed.
  ##
  ## What is true: on this route the declared thetas reach the model only
  ## through the `rxEdA.*` anchors, which nothing reads, so the observation
  ## likelihood does not depend on them.  What was wrong was concluding the
  ## family M-step is therefore mandatory.  The objective that DOES identify
  ## them is the eta density, and saem has had it all along -- it was trapped
  ## inside `etaDistMstep()`, whose loop is gated on that control.
  ##
  ## So the control now selects nothing here: the family MLE fits native
  ## parameters and inverts them, while Q2 maximizes the same likelihood over
  ## the thetas directly.  Letting both run measured MARE 17.98% against Q2's
  ## 3.27%, so `etaDistMstep` is forced inert on this route and the two settings
  ## must agree EXACTLY.
  skip_on_cran()
  .d <- .edDirectData(n = 40)
  .est <- function(.ms) {
    .f <- suppressWarnings(nlmixr2(
      .edDirectModel(), .d, "saem",
      saemControl(nBurn = 15, nEm = 15, nmc = 3, print = 0, seed = 99,
                  etaDistParam = "direct", etaDistMstep = .ms)))
    setNames(.f$parFixedDf$Estimate, rownames(.f$parFixedDf))
  }
  .on <- .est(TRUE)
  .off <- .est(FALSE)
  for (.n in c("lclm", "lclrv", "lv", "prop.sd")) {
    expect_equal(unname(.on[[.n]]), unname(.off[[.n]]), tolerance = 1e-10,
                 info = .n)
  }
  ## and neither is frozen at ini()
  expect_gt(abs(.on[["lclm"]] - log(4)), 0.05)
  ## the cdf route is unaffected -- there etaDistMstep still selects a real
  ## alternative, and the two settings must NOT agree
  expect_error(.preProcessEtaDist(rxode2::as.rxUi(.edDirectModel()), est = "saem",
                                  control = saemControl(etaDistParam = "cdf",
                                                        etaDistMstep = FALSE,
                                                        etaDistWarmStart = FALSE)),
               NA)
})

test_that("rxEtaDistExpand refuses what the direct prior cannot express", {
  ## a declared eta correlated with an ORDINARY one: the prior splits into a
  ## family part and a Gaussian part, and that split is exact only when the two
  ## do not share an omega block
  .f <- function() {
    ini({
      lclm <- log(4)
      lv <- log(50)
      lclrv <- log(0.8)
      prop.sd <- c(0, 0.2)
      dist(eta.cl) ~ dgamma(shape = 1/exp(lclrv),
                            rate = 1/(exp(lclrv)*exp(lclm)))
      eta.v ~ 0.1
    })
    model({
      cl <- eta.cl
      v <- exp(lv + eta.v)
      d/dt(centr) <- -cl/v*centr
      cp <- centr/v
      cp ~ prop(prop.sd)
    })
  }
  ## uncorrelated is fine
  expect_error(rxode2::rxEtaDistExpand(rxode2::as.rxUi(.f), param = "direct"), NA)
  ## correlated with an ordinary eta is not
  .g <- function() {
    ini({
      lclm <- log(4)
      lv <- log(50)
      lclrv <- log(0.8)
      prop.sd <- c(0, 0.2)
      dist(eta.cl) ~ dgamma(shape = 1/exp(lclrv),
                            rate = 1/(exp(lclrv)*exp(lclm)))
      eta.cl + eta.v ~ c(1, 0.05, 0.1)
    })
    model({
      cl <- eta.cl
      v <- exp(lv + eta.v)
      d/dt(centr) <- -cl/v*centr
      cp <- centr/v
      cp ~ prop(prop.sd)
    })
  }
  expect_error(rxode2::rxEtaDistExpand(rxode2::as.rxUi(.g), param = "direct"),
               "cannot correlate the declared")
  ## and the cdf route takes it, since there both really are normal latents
  expect_error(rxode2::rxEtaDistExpand(rxode2::as.rxUi(.g), param = "cdf"), NA)
})

test_that("the Q1/Q2 partition sends each route's thetas to the right owner", {
  ## NoLimits.jl's rule (_partition_q1_q2_names, src/estimation/common.jl:4557)
  ## applied to our model text: a theta is Q2 when it appears in a random-effect
  ## distribution expression and in NO observation-side one.
  ##
  ## The two routes land on OPPOSITE sides of it, which is the whole reason one
  ## rule serves both:
  ##
  ##   cdf     the decoder reads the rxEdA.* anchors  -> thetas are in the
  ##           observation path                       -> Q1
  ##   direct  nothing reads them                     -> prior-only  -> Q2
  ##
  ## An empty Q2 set on the cdf route is the assertion that matters most.  There
  ## the latent is a fixed N(0,1), so log p(z) is theta-free and the eta-density
  ## objective has a fixed point at the current theta -- measured, using it on a
  ## cdf model costs MARE 2.34% -> 17.57%.
  .u0 <- rxode2::as.rxUi(.edDirectModel())
  .all <- c("lclm", "lv", "lclrv", "prop.sd")
  for (.p in c("cdf", "direct")) {
    .ui <- rxode2::rxUiDecompress(
      .preProcessEtaDist(.u0, est = "saem",
                         control = saemControl(etaDistParam = .p,
                                               etaDistMstep = TRUE,
                                               etaDistWarmStart = FALSE))$ui)
    .s <- .etaDistThetaSplit(.ui, .all)
    if (.p == "cdf") {
      expect_equal(.s$q2, character(0))
      expect_setequal(.s$q1, .all)
    } else {
      expect_setequal(.s$q2, c("lclm", "lclrv"))
      expect_setequal(.s$q1, c("lv", "prop.sd"))
    }
    ## and it reaches the estimator's metadata, not just the helper
    .en <- .ui$iniDf$name[!is.na(.ui$iniDf$neta1) &
                            .ui$iniDf$neta1 == .ui$iniDf$neta2]
    .edi <- .etaDistMstepInfo(.ui, .ui$saemEtaTrans, .en, .ui$saemParamsToEstimate)
    expect_false(is.null(.edi))
    expect_equal(sum(unlist(.edi$q2)), if (.p == "direct") 2L else 0L)
  }
})

test_that("a prior-only theta is estimated with the family M-step OFF", {
  ## The engagement test for Q2, and the one that fails if it is inert.
  ##
  ## On the direct route the declared thetas reach the model only through the
  ## rxEdA.* anchors, which nothing reads, so the observation likelihood does
  ## not depend on them.  Before Q2 was reachable independently of
  ## etaDistMstep, this fit returned lclm at 1.3901 from a start of
  ## log(4) = 1.3863 -- frozen -- while prop.sd landed on 0.1518 against a truth
  ## of 0.15 and the eta sample matched the true gamma.  Nothing in that output
  ## said the declaration had not been estimated, which is why the assertion is
  ## on MOVEMENT and not only on closeness.
  skip_on_cran()
  .fit <- suppressWarnings(nlmixr2(
    .edDirectModel(), .edDirectData(n = 60), "saem",
    saemControl(nBurn = 60, nEm = 60, nmc = 3, print = 0, seed = 99,
                etaDistParam = "direct", etaDistMstep = FALSE)))
  .p <- setNames(.fit$parFixedDf$Estimate, rownames(.fit$parFixedDf))
  ## it MOVED off ini() -- 0.004 was the frozen signature
  expect_gt(abs(.p[["lclm"]] - log(4)), 0.05)
  expect_gt(abs(.p[["lclrv"]] - log(0.8)), 0.05)
  ## and it moved toward the truth, not merely away from the start
  expect_lt(abs(.p[["lclm"]] - log(5.104)), 0.25)
})
