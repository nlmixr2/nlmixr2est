# Phase 3.2: a covariate on a declaration argument used to disable the
# declared-distribution M-step for the WHOLE model.
#
# The refusal was implicit and total.  .etaDistMstepCore() evaluates each
# declaration's arguments in an environment holding only the thetas, so an
# argument reading a covariate column evaluated to NA, and
# `if (any(vapply(.args, anyNA, logical(1)))) return(NULL)` then returned no
# metadata at all.  A model with one covariate-carrying declaration and three
# plain ones lost the M-step on all four -- and silently, because NULL metadata
# downstream is indistinguishable from "this model has no declaration".
#
# The classification is now per declaration, so the plain ones are still fitted
# and only the covariate-carrying one stands down.  Its thetas then have to stay
# in the OUTER optimizer: a theta held out with nothing left to update it sits
# at its ini() value for the whole fit and is reported as an estimate.

.edcModel <- function(clRate, v1Decl) {
  ## Declare ONLY the thetas the two declarations actually reference: rxode2
  ## rejects a model with an ini() parameter the model block never uses, and
  ## these fixtures deliberately vary which thetas appear (the shared-theta case
  ## drops lv1rv, the no-covariate case drops bWT).
  .opt <- c("lclm", "lv1m", "lclrv", "lv1rv", "bWT")
  .txt <- paste(clRate, v1Decl)
  .use <- .opt[vapply(.opt, function(.z) grepl(.z, .txt, fixed = TRUE), logical(1))]
  .val <- c(lclm = "1.5", lv1m = "1.5", lclrv = "-1", lv1rv = "-1", bWT = "0.1")
  .ini <- paste0("      ", .use, " <- ", .val[.use], collapse = "\n")
  eval(parse(text = paste0("function() {
    ini({
      lka <- 0.5
", .ini, "
      eta.cl + eta.v ~ c(1, 0.3, 1)
      dist(eta.cl) ~ dgamma(shape = 1/exp(lclrv), rate = ", clRate, ")
      ", v1Decl, "
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
  }")))
}

## the ui the M-step actually sees: expanded, carrying the declaration stash
## that the expansion would otherwise destroy (same idiom as
## test-etaDistThetaSens.R -- .etaDistMstepInfoFocei() needs the rxz.* latents
## that only expansion creates, and the stash that only the hook preserves)
.edcUi <- function(f) {
  .ui <- rxode2::rxUiDecompress(nlmixr2est::nlmixr2(f))
  .st <- nlmixr2est:::.etaDistDeclStash(.ui, rxode2::rxUiEtaDists(.ui))
  .u2 <- rxode2::rxUiDecompress(rxode2::rxEtaDistExpand(.ui))
  nlmixr2est:::.etaDistDeclSet(.u2, .st)
  .u2
}

.edcCore <- function(f) nlmixr2est:::.etaDistMstepCore(.edcUi(f))
.edcInfo <- function(f) nlmixr2est:::.etaDistMstepInfoFocei(.edcUi(f))

.edcPlainCl <- "1/(exp(lclrv)*exp(lclm))"
.edcCovCl   <- "1/(exp(lclrv)*exp(lclm + bWT*(WT - 70)))"
.edcPlainV1 <- "dist(eta.v) ~ dgamma(shape = 1/exp(lv1rv), rate = 1/(exp(lv1rv)*exp(lv1m)))"
.edcCovV1   <- "dist(eta.v) ~ dgamma(shape = 1/exp(lv1rv), rate = 1/(exp(lv1rv)*exp(lv1m + bWT*(WT - 70))))"

test_that("no covariate: every declaration is usable (regression guard T1)", {
  .c <- .edcCore(.edcModel(.edcPlainCl, .edcPlainV1))
  expect_false(is.null(.c))
  expect_identical(.c$hasCov, c(FALSE, FALSE))
  expect_identical(.c$usable, c(TRUE, TRUE))
})

test_that("a covariate on ONE declaration leaves the OTHER fittable", {
  .c <- .edcCore(.edcModel(.edcCovCl, .edcPlainV1))
  # this is the whole point: metadata at all, where before it was NULL
  expect_false(is.null(.c))
  expect_identical(.c$hasCov, c(TRUE, FALSE))
  expect_identical(.c$usable, c(FALSE, TRUE))
})

test_that("a covariate declaration's thetas stay in the outer optimizer", {
  .i <- .edcInfo(.edcModel(.edcCovCl, .edcPlainV1))
  skip_if(is.null(.i), "focei metadata unavailable on this model")
  # cl's declaration carries the covariate, so the M-step does not own its
  # thetas -- holding them out would freeze them at ini() with nothing to move
  # them, which reports as an estimate and is exactly the silent failure here
  expect_false("lclm" %in% .i$thetaNames)
  expect_false("bWT" %in% .i$thetaNames)
  # v1's declaration is plain, so the M-step does own those
  expect_true("lv1m" %in% .i$thetaNames)
  expect_true("lv1rv" %in% .i$thetaNames)
  expect_identical(as.integer(.i$usable), c(0L, 1L))
})

test_that("a theta SHARED with a covariate declaration is not held out", {
  # both declarations read lclrv; only v1's is fittable.  Holding lclrv out for
  # v1's sake would freeze it for cl's covariate declaration too.
  .shared <- "dist(eta.v) ~ dgamma(shape = 1/exp(lclrv), rate = 1/(exp(lclrv)*exp(lv1m)))"
  .i <- .edcInfo(.edcModel(.edcCovCl, .shared))
  skip_if(is.null(.i), "focei metadata unavailable on this model")
  expect_false("lclrv" %in% .i$thetaNames)
  expect_true("lv1m" %in% .i$thetaNames)
})

test_that("the map closure declines LOUDLY for a declaration it does not own", {
  .i <- .edcInfo(.edcModel(.edcCovCl, .edcPlainV1))
  skip_if(is.null(.i), "focei metadata unavailable on this model")
  # there is no single population `a` to invert when an argument varies by
  # subject, so a best-effort answer would be a wrong number, not a missing one
  expect_null(.i$map(1L, c(1, 1)))
  expect_false(is.null(.i$map(2L, as.numeric(.i$args[2, seq_len(2)]))))
})

test_that("ALL declarations covariate-carrying still stands the M-step down", {
  .c <- .edcCore(.edcModel(.edcCovCl, .edcCovV1))
  expect_null(.c)
})

test_that("the warm start NAMES the covariate instead of refusing silently", {
  .ui <- rxode2::rxUiDecompress(nlmixr2est::nlmixr2(.edcModel(.edcCovCl, .edcPlainV1)))
  expect_warning(.s <- nlmixr2est:::.etaDistSurrogate(.ui), "WT")
  expect_null(.s)
})

test_that("covRef lets the warm start proceed", {
  .ui <- rxode2::rxUiDecompress(nlmixr2est::nlmixr2(.edcModel(.edcCovCl, .edcPlainV1)))
  expect_silent(.s <- nlmixr2est:::.etaDistSurrogate(.ui, covRef = list(WT = 70)))
  expect_false(is.null(.s))
})

test_that("covRef evaluates the moments AT the reference", {
  # at WT = 70 the covariate term vanishes, so the moments must equal the
  # plain declaration's -- a reference that was ignored would not
  .cov <- "dgamma(shape = 1/exp(lclrv), rate = 1/(exp(lclrv)*exp(lclm + bWT*(WT - 70))))"
  .plain <- "dgamma(shape = 1/exp(lclrv), rate = 1/(exp(lclrv)*exp(lclm)))"
  .tv <- list(lclrv = -2.0, lclm = 1.5, bWT = 0.1)
  .a <- nlmixr2est:::.etaDistMoments(.cov, .tv, covRef = list(WT = 70))
  .b <- nlmixr2est:::.etaDistMoments(.plain, .tv)
  expect_equal(.a[["mean"]], .b[["mean"]], tolerance = 1e-10)
  expect_equal(.a[["var"]], .b[["var"]], tolerance = 1e-10)
  # and a DIFFERENT reference must move it, or the argument is being dropped
  .c <- nlmixr2est:::.etaDistMoments(.cov, .tv, covRef = list(WT = 90))
  expect_false(isTRUE(all.equal(.a[["mean"]], .c[["mean"]])))
})

test_that("a theta is not shadowed by a covariate of the same name", {
  .tv <- list(lclrv = -2.0, lclm = 1.5)
  .plain <- "dgamma(shape = 1/exp(lclrv), rate = 1/(exp(lclrv)*exp(lclm)))"
  .a <- nlmixr2est:::.etaDistMoments(.plain, .tv)
  .b <- nlmixr2est:::.etaDistMoments(.plain, .tv, covRef = list(lclm = 99))
  expect_equal(.a[["mean"]], .b[["mean"]], tolerance = 1e-12)
})

# --------------------------------------------------------- end to end (T3) --
# The metadata tests above work on a hand-built ui.  This one goes through the
# real path -- nlmixr2() then rxEtaDistExpand() -- and checks the thing that
# actually makes a covariate on a declaration work at all: the expansion puts
# the covariate INTO the decoder line, so the ordinary per-record solve
# evaluates it and no M-step machinery is needed for it.

.edcCovModel <- function() {
  function() {
    ini({ lclm <- 1.5; lv1m <- 1.5; lclrv <- -1.2; lv1rv <- -1.2; bWT <- 0.5
          prop.sd <- 0.1
          eta.cl + eta.v ~ c(1, 0.3, 1) })
    model({
      dist(eta.cl) ~ dgamma(shape = 1/exp(lclrv),
                            rate = 1/(exp(lclrv)*exp(lclm + bWT*log(WT/70))))
      dist(eta.v)  ~ dgamma(shape = 1/exp(lv1rv), rate = 1/(exp(lv1rv)*exp(lv1m)))
      cl <- eta.cl; v <- eta.v
      d/dt(central) <- -(cl/v)*central
      cp <- central/v
      cp ~ prop(prop.sd)
    })
  }
}

test_that("the covariate lands in the decoder line, per record", {
  skip_on_cran()
  .ui <- rxode2::rxUiDecompress(nlmixr2est::nlmixr2(.edcCovModel()))
  expect_true("WT" %in% .ui$allCovs)
  .e <- rxode2::rxUiDecompress(rxode2::rxEtaDistExpand(.ui))
  .ln <- vapply(.e$lstExpr, function(.z) paste(deparse(.z), collapse = " "),
                character(1))
  .dec <- .ln[grepl("^eta.cl <-", .ln)]
  expect_length(.dec, 1L)
  # this is what makes a covariate on a declaration need NO M-step machinery:
  # the solve evaluates the argument per record, covariate and all
  expect_true(grepl("WT", .dec, fixed = TRUE))
  expect_true(grepl("gammapInv", .dec, fixed = TRUE))
  # and the declaration WITHOUT the covariate does not acquire one
  expect_false(grepl("WT", .ln[grepl("^eta.v <-", .ln)], fixed = TRUE))
})

test_that("through the real path, only the plain declaration is held out", {
  skip_on_cran()
  .ui <- rxode2::rxUiDecompress(nlmixr2est::nlmixr2(.edcCovModel()))
  .st <- nlmixr2est:::.etaDistDeclStash(.ui, rxode2::rxUiEtaDists(.ui))
  .e <- rxode2::rxUiDecompress(rxode2::rxEtaDistExpand(.ui))
  nlmixr2est:::.etaDistDeclSet(.e, .st)
  .c <- nlmixr2est:::.etaDistMstepCore(.e)
  expect_false(is.null(.c))
  expect_identical(.c$hasCov, c(TRUE, FALSE))
  .i <- nlmixr2est:::.etaDistMstepInfoFocei(.e)
  expect_false(is.null(.i))
  # the covariate declaration's thetas stay estimable by the outer optimizer
  for (.t in c("lclm", "lclrv", "bWT")) expect_false(.t %in% .i$thetaNames)
  # the plain one's, and the copula, are owned by the M-step
  for (.t in c("lv1m", "lv1rv")) expect_true(.t %in% .i$thetaNames)
  expect_true(any(grepl("^rxCor[.]", .i$thetaNames)))
})

test_that("the M-step warns rather than silently mis-estimating a covariate", {
  skip_on_cran()
  # MEASURED: with the M-step ON, focei returned bWT -0.027 for a truth of
  # +0.75; with it OFF, +0.598, against +0.574 from a log-normal reference fit
  # of the same data.  The cause is upstream of the M-step -- the C++ argument
  # parser resolves symbols against the declaration's THETA names only, so an
  # argument reading a data column fails to parse (etaDistExprParse,
  # src/etaDistExpr.h) and the general family objective never runs.  Until the
  # covariate names are passed to it (nSym/rec, already accepted by
  # rxEtaDistLoglikObj), this route must SAY so rather than return a number.
  .ui <- rxode2::rxUiDecompress(nlmixr2est::nlmixr2(.edcModel(.edcCovCl, .edcPlainV1)))
  .st <- nlmixr2est:::.etaDistDeclStash(.ui, rxode2::rxUiEtaDists(.ui))
  .u2 <- rxode2::rxUiDecompress(rxode2::rxEtaDistExpand(.ui))
  nlmixr2est:::.etaDistDeclSet(.u2, .st)
  rxode2::rxAssignControlValue(.u2, "etaDistMstep", TRUE)
  expect_warning(nlmixr2est:::.foceiEtaDistSetup(.u2), "covariate coefficient")
  expect_warning(nlmixr2est:::.foceiEtaDistSetup(.u2), "eta.cl")
  expect_warning(nlmixr2est:::.foceiEtaDistSetup(.u2), "etaDistMstep=FALSE")
})

test_that("no covariate means no warning", {
  skip_on_cran()
  .ui <- rxode2::rxUiDecompress(nlmixr2est::nlmixr2(.edcModel(.edcPlainCl, .edcPlainV1)))
  .st <- nlmixr2est:::.etaDistDeclStash(.ui, rxode2::rxUiEtaDists(.ui))
  .u2 <- rxode2::rxUiDecompress(rxode2::rxEtaDistExpand(.ui))
  nlmixr2est:::.etaDistDeclSet(.u2, .st)
  rxode2::rxAssignControlValue(.u2, "etaDistMstep", TRUE)
  # every declaration usable -> the guard must stay quiet
  expect_no_warning(nlmixr2est:::.foceiEtaDistSetup(.u2))
})

test_that("the argument parser is what cannot see the covariate", {
  # pins the exact mechanism, so a fix is visible here first: the grammar
  # already handles the expression -- only the SYMBOL is unknown
  .f <- nlmixr2est:::rxEtaDistArgsToThetasTest_
  .plain <- .f(c("1/exp(lclrv)", "1/(exp(lclrv)*exp(lclm))"),
               c("lclrv", "lclm"), c(-2, 1.5), c(11, 2.2))
  expect_false(is.null(.plain))
  expect_length(.plain, 2L)
  # same shape of expression, plus a covariate symbol -> declines
  .cov <- .f(c("1/exp(lclrv)", "1/(exp(lclrv)*exp(lclm + bWT*log(WT/70)))"),
             c("lclrv", "lclm", "bWT"), c(-2, 1.5, 0.2), c(11, 2.2))
  expect_true(is.null(.cov) || length(.cov) == 0L)
  # and it parses as soon as the symbol is known, which is the whole fix
  .named <- .f(c("1/exp(lclrv)", "1/(exp(lclrv)*exp(lclm + bWT*log(WT/70)))"),
               c("lclrv", "lclm", "bWT", "WT"), c(-2, 1.5, 0.2, 70), c(11, 2.2))
  expect_false(is.null(.named))
  expect_length(.named, 4L)
})

# ---------------------------- the EBE route estimates a covariate directly --
# Each subject's EBE for a declared parameter is a draw from family(args_i)
# with args_i depending on THAT subject's covariates, so the family likelihood
# over the EBE sample identifies the covariate coefficient -- no
# observation-likelihood route needed.  rxEtaDistLoglikObj() implements exactly
# that: it maximizes sum_r wt[r]*log p(eta[r]; args_r(theta, rec_r)) over the
# THETAS, instead of fitting one population native parameter set and inverting
# it (there is no single population `a` to invert when an argument varies by
# subject).  These pin the estimator itself, away from any estimator's gating.

test_that("the per-record objective recovers a covariate from the EBEs", {
  skip_on_cran()
  set.seed(7)
  .f <- nlmixr2est:::rxEtaDistLoglikTest_
  .GAMMA <- 13L
  .n <- 4000L; .lclrv <- -2.4; .lclm <- 1.63; .bWT <- 0.75
  .WT <- stats::runif(.n, 40, 100)
  # drawn from the model the objective assumes, so nothing but the covariate
  # can explain the effect
  .eta <- stats::rgamma(.n, shape = 1/exp(.lclrv),
                        rate = 1/(exp(.lclrv)*exp(.lclm + .bWT*log(.WT/70))))
  .exprs <- c("1/exp(lclrv)", "1/(exp(lclrv)*exp(lclm + bWT*log(WT/70)))")
  .vars <- c("lclrv", "lclm", "bWT", "WT")   # thetas THEN the record symbols
  .rec <- matrix(.WT, ncol = 1L); .wt <- rep(1, .n)
  .obj <- function(.th) {
    .v <- .f(.GAMMA, .exprs, .vars, .th, .rec, .eta, .wt)
    if (length(.v) == 0L) NA_real_ else -.v[1]
  }
  # the truth is a better explanation than no covariate, by a wide margin
  expect_lt(.obj(c(.lclrv, .lclm, .bWT)), .obj(c(.lclrv, .lclm, 0)))
  .o <- stats::optim(c(-2.0, 1.5, 0.2), .obj, method = "Nelder-Mead",
                     control = list(maxit = 4000, reltol = 1e-10))
  expect_equal(.o$par[1], .lclrv, tolerance = 0.05)
  expect_equal(.o$par[2], .lclm, tolerance = 0.05)
  expect_equal(.o$par[3], .bWT, tolerance = 0.08)
})

test_that("the objective declines when the covariate symbol is not supplied", {
  # same contract as the parser: an unresolvable symbol falls back rather than
  # guessing, which is what silently cost the covariate before
  .f <- nlmixr2est:::rxEtaDistLoglikTest_
  .v <- .f(13L, c("1/exp(lclrv)", "1/(exp(lclrv)*exp(lclm + bWT*log(WT/70)))"),
           c("lclrv", "lclm", "bWT"), c(-2.4, 1.63, 0.75),
           matrix(0, nrow = 3L, ncol = 0L), c(1, 2, 3), c(1, 1, 1))
  expect_length(.v, 0L)
})

test_that("nSym == 0 still works, so no-covariate models are unchanged", {
  .f <- nlmixr2est:::rxEtaDistLoglikTest_
  .eta <- c(2.0, 3.0, 5.0, 7.0)
  .v <- .f(13L, c("1/exp(lclrv)", "1/(exp(lclrv)*exp(lclm))"),
           c("lclrv", "lclm"), c(-2.4, 1.63),
           matrix(0, nrow = 4L, ncol = 0L), .eta, rep(1, 4))
  expect_length(.v, 1L)
  # it is the plain gamma log-likelihood, so check it against dgamma directly
  .sh <- 1/exp(-2.4); .rt <- 1/(exp(-2.4)*exp(1.63))
  expect_equal(.v[1], sum(stats::dgamma(.eta, shape = .sh, rate = .rt, log = TRUE)),
               tolerance = 1e-8)
})

# ------------------------- the phi mapping must not count covariate coefficients --
# saemParamsToEstimate INTERLEAVES each theta with its mu-referenced covariate
# coefficients, because it indexes MCOV.  A coefficient is not a phi parameter
# (covstruct stays 7x7 while that list grows to 8), but the indices built from
# it are consumed as PHI indices.  Matching against the longer list put every
# parameter at or after the coefficient one too high, so the M-step wrote into
# the WRONG slots: declaration 1's (lclrv, lclm, bWT) landed on phi0 columns
# (3, 0, 1) = (lv1rv, lclm, lv1m).  A range check cannot catch it -- the indices
# stay within nphi, they just mean something else -- so pin the mapping itself.

test_that("no declared theta comes back holding another parameter's start value", {
  skip_on_cran()
  skip_on_os("windows")
  # THE regression, and it has to go through a fit: saemParamsToEstimate only
  # interleaves the covariate coefficient once saem's own setup has classified
  # it, so neither a fresh nor an expanded ui reproduces the condition -- both
  # report an EMPTY saemMuRefCovariateDataFrame.
  #
  # With the bug the M-step wrote into the wrong phi0 columns and the table came
  # back with lclrv holding lv1m's start value and the copula theta holding
  # lv1rv's -- EXACT equalities, which is what makes this testable.  Distinct
  # starting values are the whole trick.
  .d <- suppressWarnings(nlmixr2data::theo_sd)
  skip_if(is.null(.d), "theo_sd unavailable")
  .d$WT <- 60 + (.d$ID %% 5L) * 8      # subject-constant, so saem accepts it
  .m <- function() {
    ini({ lclm <- 1.11; lv1m <- 2.22; lclrv <- -3.33; lv1rv <- -4.44; bWT <- 0.55
          prop.sd <- 0.77
          eta.cl + eta.v ~ c(1, 0.3, 1) })
    model({
      dist(eta.cl) ~ dgamma(shape = 1/exp(lclrv),
                            rate = 1/(exp(lclrv)*exp(lclm + bWT*log(WT/70))))
      dist(eta.v)  ~ dgamma(shape = 1/exp(lv1rv), rate = 1/(exp(lv1rv)*exp(lv1m)))
      cl <- eta.cl; v <- eta.v
      d/dt(central) <- -(cl/v)*central
      cp <- central/v
      cp ~ prop(prop.sd)
    })
  }
  .f <- try(suppressWarnings(suppressMessages(nlmixr2est::nlmixr2(
    .m, .d, est = "saem",
    control = nlmixr2est::saemControl(print = 0, covMethod = "", calcTables = FALSE,
                                      etaDistMstep = TRUE, nBurn = 40, nEm = 40,
                                      seed = 99)))), silent = TRUE)
  skip_if(inherits(.f, "try-error"), "saem fit unavailable in this environment")
  .e <- .f$parFixedDf[, "Estimate"]
  names(.e) <- rownames(.f$parFixedDf)
  skip_if(!all(c("lclrv", "rxCor.eta.v.eta.cl") %in% names(.e)), "parameters absent")
  # the two the bug corrupted, against the values it handed them
  expect_false(isTRUE(all.equal(unname(.e[["lclrv"]]), 2.22)))
  expect_false(isTRUE(all.equal(unname(.e[["rxCor.eta.v.eta.cl"]]), -4.44)))
})

test_that("focei's M-step honors etaDistEvery and fires once per outer evaluation", {
  skip_on_cran()
  expect_equal(nlmixr2est::foceiControl()$etaDistEvery, 20L)
  expect_equal(nlmixr2est::foceiControl(etaDistEvery = 5L)$etaDistEvery, 5L)
  # saem's cadence and focei's should agree by default -- they are the same knob
  expect_equal(nlmixr2est::foceiControl()$etaDistEvery,
               nlmixr2est::saemControl()$etaDistEvery)
})

test_that("a coarser cadence fires the M-step far less often", {
  skip_on_cran()
  skip_on_os("windows")
  # The step rewrites the thetas it owns, so firing it inside the objective made
  # the objective at the same theta differ between calls and every outer
  # optimizer stalled.  It used to run once per pass of the inner
  # {re-optimize etas, update} loop; it now runs at most once per OUTER
  # evaluation, thinned by etaDistEvery.  Counting the firings is the direct
  # observable -- ~200 vs a handful on the same fit.
  .d <- suppressWarnings(nlmixr2data::theo_sd)
  skip_if(is.null(.d), "theo_sd unavailable")
  .d$WT <- 60 + (.d$ID %% 5L) * 8
  .m <- function() {
    ini({ lclm <- 1.5; lv1m <- 1.5; lclrv <- -1.2; lv1rv <- -1.2; bWT <- 0.2
          prop.sd <- 0.15
          eta.cl + eta.v ~ c(1, 0.3, 1) })
    model({
      dist(eta.cl) ~ dgamma(shape = 1/exp(lclrv),
                            rate = 1/(exp(lclrv)*exp(lclm + bWT*log(WT/70))))
      dist(eta.v)  ~ dgamma(shape = 1/exp(lv1rv), rate = 1/(exp(lv1rv)*exp(lv1m)))
      cl <- eta.cl; v <- eta.v
      d/dt(central) <- -(cl/v)*central
      cp <- central/v
      cp ~ prop(prop.sd)
    })
  }
  .fire <- function(.every) {
    .f <- try(suppressWarnings(suppressMessages(nlmixr2est::nlmixr2(
      .m, .d, est = "focei",
      control = nlmixr2est::foceiControl(print = 0, covMethod = "",
                                         calcTables = FALSE, etaDistMstep = TRUE,
                                         etaDistEvery = .every)))), silent = TRUE)
    if (inherits(.f, "try-error")) return(NA_integer_)
    tryCatch(nlmixr2est:::foceiEtaDistN_(), error = function(e) NA_integer_)
  }
  .n1 <- .fire(1L)
  .n20 <- .fire(20L)
  skip_if(is.na(.n1) || is.na(.n20), "focei fit unavailable in this environment")
  # the cadence has to actually thin it, not merely be accepted
  expect_gt(.n1, .n20)
})

test_that("the support rejection is tunable on both controls", {
  # Threshold: 0 rejects only an exact zero; raising it also rejects a merely
  # SUBNORMAL draw, which is not a draw either.
  expect_equal(nlmixr2est::foceiControl()$etaDistSupportEps, 0)
  expect_equal(nlmixr2est::impmapControl()$etaDistSupportEps, 0)
  expect_equal(nlmixr2est::impmapControl(etaDistSupportEps = 1e-300)$etaDistSupportEps,
               1e-300)
  expect_equal(nlmixr2est::foceiControl(etaDistSupportEps = 1e-300)$etaDistSupportEps,
               1e-300)
  # a negative threshold would reject nothing AND admit negative etas into a
  # positive-support family, so refuse it rather than silently disabling
  expect_error(nlmixr2est::impmapControl(etaDistSupportEps = -1))
  # and it round-trips through do.call(), which is how the control is rebuilt
  .c <- nlmixr2est::impmapControl(etaDistSupportEps = 1e-300)
  expect_equal(do.call(nlmixr2est::impmapControl,
                       .c[names(.c) %in% names(formals(nlmixr2est::impmapControl))]
                       )$etaDistSupportEps, 1e-300)
})
