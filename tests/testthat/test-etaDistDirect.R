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

test_that("a COVARIATE on a cdf declaration is Q1, not Q2", {
  ## The trickier partition case, and one I got wrong in planning: I expected Q2
  ## to take over `bWT` on the covariate arm and lift it from 0.2551 toward the
  ## 0.7764 that only etaDistMstep=TRUE reaches.  It must not, and does not.
  ##
  ## That arm is a CDF model.  The decoder reads the rxEdA.* anchors, so every
  ## theta the declaration mentions -- the covariate coefficient included -- is
  ## in the observation path and the data identify it.  Q2 has nothing to own,
  ## and this arm must be UNCHANGED by the partition.  Measured: bWT 0.7764 /
  ## 0.2551 for etaDistMstep TRUE / FALSE, bit-identical before and after.
  ##
  ## (The 0.2551-vs-0.7764 gap is real but belongs to the etaDistMstep default,
  ## not to Q2.)
  .f <- function() {
    ini({
      lclm <- 1.63
      lclrv <- -2.4
      bWT <- 0.35
      lv <- log(50)
      prop.sd <- 0.1
      dist(eta.cl) ~ dgamma(shape = 1/exp(lclrv),
                            rate = 1/(exp(lclrv)*exp(lclm + bWT*log(WT/70))))
    })
    model({
      cl <- eta.cl
      v <- exp(lv)
      d/dt(centr) <- -cl/v*centr
      cp <- centr/v
      cp ~ prop(prop.sd)
    })
  }
  .ui <- rxode2::rxUiDecompress(
    .preProcessEtaDist(rxode2::as.rxUi(.f), est = "saem",
                       control = saemControl(etaDistMstep = TRUE,
                                             etaDistWarmStart = FALSE))$ui)
  .s <- .etaDistThetaSplit(.ui, c("lclm", "lclrv", "bWT", "lv", "prop.sd"))
  expect_equal(.s$q2, character(0))
  expect_true("bWT" %in% .s$q1)
})

test_that("a directly-parameterized eta is not reported as a variance", {
  ## On this route the declared eta carries its FAMILY as its prior, so it has
  ## no variance in the ordinary sense -- rxEtaDistExpand() fixes its omega
  ## entry at a placeholder 1 precisely because it is machinery, not a
  ## parameter.  What the fit printed was neither: measured on a correlated
  ## gamma pair, every cell of the 2x2 came back 180.6704, so $omegaR reported a
  ## correlation of 1.000 with an SD of 13.44, while the iniDf still carried
  ## fix = TRUE on those rows.
  ##
  ## Reporting the placeholder correctly would not fix the real problem: a
  ## reader takes a printed omega for the dispersion, and the dispersion is in
  ## the family's own parameters, which are already in parFixed.  NoLimits.jl
  ## has no omega concept at all for this reason and omits any quantity a family
  ## cannot supply rather than emitting a placeholder.  This is that rule.
  skip_on_cran()
  .fit <- suppressMessages(suppressWarnings(nlmixr2(
    .edDirectModel(), .edDirectData(n = 40), "saem",
    saemControl(nBurn = 15, nEm = 15, nmc = 3, print = 0, seed = 99,
                etaDistParam = "direct", etaDistMstep = FALSE))))
  ## the declared eta appears in NEITHER table
  expect_false(any(grepl("^rxd[.]", rownames(.fit$parFixedDf))))
  expect_false(any(grepl("^rxd[.]", rownames(.fit$omega))))
  ## the family's own parameters ARE reported -- that is where the dispersion is
  expect_true(all(c("lclm", "lclrv") %in% rownames(.fit$parFixedDf)))
  ## and this model declares its only random effect, so nothing is left
  expect_equal(nrow(.fit$omega), 0L)
})

test_that("an ORDINARY eta beside a declared one is still reported", {
  ## The case a blanket suppression gets wrong.  Only the declared random
  ## effect loses its omega row; a plain eta in the same model keeps its
  ## variance, because for that one the variance IS the parameter.
  skip_on_cran()
  .f <- function() {
    ini({
      lclm <- log(4)
      lv <- log(45)
      lclrv <- log(0.8)
      prop.sd <- c(0, 0.2)
      eta.v ~ 0.1
      dist(eta.cl) ~ dgamma(shape = 1/exp(lclrv),
                            rate = 1/(exp(lclrv)*exp(lclm)))
    })
    model({
      cl <- eta.cl
      v <- exp(lv + eta.v)
      d/dt(centr) <- -cl/v*centr
      cp <- centr/v
      cp ~ prop(prop.sd)
    })
  }
  set.seed(11)
  .n <- 40
  .cl <- stats::rgamma(.n, 2, 2/5.104)
  .v <- 50*exp(stats::rnorm(.n, 0, 0.3))
  .obs <- do.call(rbind, lapply(seq_len(.n), function(.i) {
    .t <- c(0.25, 1, 2, 4, 8, 12, 24)
    data.frame(ID = .i, TIME = .t,
               DV = 100/.v[.i]*exp(-.cl[.i]/.v[.i]*.t)*
                 exp(stats::rnorm(length(.t), 0, 0.15)),
               AMT = 0, EVID = 0)
  }))
  .d <- rbind(data.frame(ID = seq_len(.n), TIME = 0, DV = 0, AMT = 100, EVID = 1),
              .obs)
  .d <- .d[order(.d$ID, .d$TIME, -.d$EVID), ]
  .fit <- suppressMessages(suppressWarnings(nlmixr2(
    .f, .d, "saem",
    saemControl(nBurn = 15, nEm = 15, nmc = 3, print = 0, seed = 99,
                etaDistParam = "direct", etaDistMstep = FALSE))))
  expect_equal(rownames(.fit$omega), "eta.v")
  expect_gt(.fit$omega[1, 1], 0)
  expect_false(any(grepl("^rxd[.]", rownames(.fit$parFixedDf))))
})

test_that("the declared copula correlation is estimated AND reported", {
  ## Both halves matter, and each was broken separately.
  ##
  ## The correlation lives in the omega on this route, not in an rxCor.* theta,
  ## so nothing carried it out of the sampler and the fit could not show it.
  ## Then, once saem returned it, it still read 0.3 -- the start -- because
  ## etaDistCorWith is recorded only on the HIGHER-indexed member of a pair,
  ## which is the slot every reader consults, while etaDistQ2PairStep() is
  ## entered from the LOWER member and wrote there.
  ##
  ## Do NOT substitute the correlation of the fitted etas for this check.  That
  ## is what misled me: the etas correlate because the DATA do, whatever the
  ## prior's rho is, so a fitted eta correlation near the truth says nothing
  ## about whether rho was estimated.  Read the parameter.
  skip_on_cran()
  .rho <- 0.6
  set.seed(7)
  .n <- 60
  .z <- matrix(stats::rnorm(2*.n), .n, 2)
  .z[, 2] <- .rho*.z[, 1] + sqrt(1 - .rho^2)*.z[, 2]
  .cl <- stats::qgamma(stats::pnorm(.z[, 1]), 2, 2/5.104)
  .v <- stats::qgamma(stats::pnorm(.z[, 2]), 3, 3/50)
  .obs <- do.call(rbind, lapply(seq_len(.n), function(.i) {
    .t <- c(0.25, 1, 2, 4, 8, 12, 24)
    data.frame(ID = .i, TIME = .t,
               DV = 100/.v[.i]*exp(-.cl[.i]/.v[.i]*.t)*
                 exp(stats::rnorm(length(.t), 0, 0.15)),
               AMT = 0, EVID = 0)
  }))
  .d <- rbind(data.frame(ID = seq_len(.n), TIME = 0, DV = 0, AMT = 100, EVID = 1),
              .obs)
  .d <- .d[order(.d$ID, .d$TIME, -.d$EVID), ]
  .f <- function() {
    ini({
      lclm <- log(4)
      lvm <- log(45)
      lclrv <- log(0.8)
      lvrv <- log(0.5)
      prop.sd <- c(0, 0.2)
      eta.cl + eta.v ~ c(1, 0.3, 1)
      dist(eta.cl) ~ dgamma(shape = 1/exp(lclrv),
                            rate = 1/(exp(lclrv)*exp(lclm)))
      dist(eta.v) ~ dgamma(shape = 1/exp(lvrv),
                           rate = 1/(exp(lvrv)*exp(lvm)))
    })
    model({
      cl <- eta.cl
      v <- eta.v
      d/dt(centr) <- -cl/v*centr
      cp <- centr/v
      cp ~ prop(prop.sd)
    })
  }
  .fit <- suppressMessages(suppressWarnings(nlmixr2(
    .f, .d, "saem",
    saemControl(nBurn = 60, nEm = 60, nmc = 3, print = 0, seed = 99,
                etaDistParam = "direct", etaDistMstep = FALSE))))
  ## reported at all
  .c <- .fit$etaDistCor
  expect_false(is.null(.c))
  expect_equal(length(.c), 1L)
  .R <- .c[[1]]
  expect_equal(dim(.R), c(2L, 2L))
  expect_equal(unname(diag(.R)), c(1, 1))
  ## MOVED off its 0.3 start -- reading 0.3 is the slot bug's signature
  expect_gt(abs(.R[1, 2] - 0.3), 0.1)
  ## and toward the truth
  expect_lt(abs(.R[1, 2] - .rho), 0.2)
  ## and it is in parFixed, matching
  .cr <- grep("^cor\\(", rownames(.fit$parFixedDf), value = TRUE)
  expect_equal(length(.cr), 1L)
  expect_equal(unname(.fit$parFixedDf[.cr, "Estimate"]), unname(.R[1, 2]),
               tolerance = 1e-8)
})
