# Phase 4.2 T5: the DEGENERATE covariate arm, end to end.
#
# Every other covariate test in test-etaDistCovariate.R is structural -- it
# asserts the parser sees the covariate, the decoder line carries it per record,
# the right thetas are held out.  None of them fits anything, so none of them
# would notice the estimator recovering a covariate effect that is not there.
#
# T5 simulates a covariate whose true coefficient is EXACTLY zero.  The obvious
# test -- fit it and assert bWT ~ 0 -- is WORTHLESS on its own, and that is the
# main thing this file records: the coefficient's start is also 0, so an
# estimator that never moves it passes.
#
# What the arm shows once the starts are displaced instead:
#
#   * focei DOES recover the coefficient -- 0.7353 against a true 0.75 -- but
#     only from a non-zero slope start and with the residual variance bounded
#     below.  From a start of exactly 0 it returns 0.0017 and stops.
#
#     Replicated at 120 subjects over independent simulations (truth 0.75):
#     with the covariate FIXED per subject the mean is 0.7132 (sd 0.061, n=5),
#     i.e. unbiased; with it varying WITHIN subject the mean is 0.6767
#     (sd 0.100, n=12), about 10% low at 2.5 standard errors of the mean.  So
#     the time-varying path carries a mild downward attenuation, and single
#     realizations scatter widely -- 0.44 and 0.80 both occur -- which is worth
#     knowing before reading any one fit as evidence of a defect.
#   * saem, on the zero-effect arm with the other thetas displaced, INVENTS an
#     effect: bWT = -1.0876, and pushes lclrv from -2.0 to -1.4139 (away from
#     the true -2.4).  On the arm where the effect is real it instead freezes
#     that declaration's thetas bit-exactly at their starts.
#
# The saem behavior is the failure this arm exists to catch -- a method that
# reports an effect where there is none selects spurious covariates on real
# data -- and it is invisible when every parameter starts at its true value,
# where saem reproduces its ini() to four decimals and looks perfect.
#
# The data is built inline from a fixed seed rather than read from
# inst/sim/simCovT45.R's output, so the test is self-contained.  It is the same
# construction at a smaller subject count: gamma CL and V1 through a Gaussian
# copula, the covariate entering the RATE of the cl declaration per record.
.edT5Data <- function(nSub = 60L, bWT = 0.0, seed = 20260912L) {
  set.seed(seed)
  .lclm <- 1.63; .lv1m <- 1.55; .lclrv <- -2.4; .lv1rv <- -2.4; .rho <- 0.5
  .tim <- c(0.25, 0.5, 1, 2, 4, 8, 12, 24)
  .u <- function(z) pmin(pmax(stats::pnorm(z), 1e-15), 1 - 1e-15)
  .z1 <- stats::rnorm(nSub); .z2 <- stats::rnorm(nSub)
  .w2 <- .rho * .z1 + sqrt(1 - .rho^2) * .z2
  .wtBase <- stats::rnorm(nSub, 70, 12)
  .wtRec <- lapply(seq_len(nSub), function(i)
    round(.wtBase[i] + cumsum(stats::rnorm(length(.tim), 0, 1.5)), 1))
  .shCL <- 1 / exp(.lclrv); .shV1 <- 1 / exp(.lv1rv)
  .V1 <- stats::qgamma(.u(.w2), shape = .shV1,
                       rate = 1 / (exp(.lv1rv) * exp(.lv1m)))
  .m <- rxode2::rxode2({
    d/dt(central) <- -(cl / v) * central
    cp <- central / v
  })
  .rows <- lapply(seq_len(nSub), function(i) {
    .cl <- stats::qgamma(.u(.z1[i]), shape = .shCL,
                         rate = 1 / (exp(.lclrv) *
                                     exp(.lclm + bWT * log(.wtRec[[i]] / 70))))
    .ev <- data.frame(id = i, time = c(0, .tim),
                      amt = c(100, rep(NA_real_, length(.tim))),
                      evid = c(1L, rep(0L, length(.tim))), cmt = 1L,
                      cl = c(.cl[1], .cl), v = .V1[i])
    .s <- rxode2::rxSolve(.m, .ev, returnType = "data.frame")
    .s <- .s[!is.na(.s$cp) & .s$time > 0, ]
    data.frame(ID = i, TIME = .s$time, CP = .s$cp, AMT = NA_real_,
               EVID = 0L, CMT = 1L,
               WT = .wtRec[[i]][seq_len(nrow(.s))])
  })
  .obs <- do.call(rbind, .rows)
  # assay limit on the TRUE concentration: without it a 1e-12 observation
  # against a 1e-5 prediction is a relative residual of 1e7 under prop(), and
  # those records dominate the residual parameter
  .obs <- .obs[.obs$CP > 0.01, ]
  .obs$DV <- .obs$CP * (1 + stats::rnorm(nrow(.obs), 0, 0.10))
  .obs <- .obs[.obs$DV > 0, c("ID", "TIME", "DV", "AMT", "EVID", "CMT", "WT")]
  .dose <- data.frame(ID = seq_len(nSub), TIME = 0, DV = NA_real_, AMT = 100,
                      EVID = 1L, CMT = 1L, WT = round(.wtBase, 1))
  .d <- rbind(.dose, .obs)
  .d[order(.d$ID, .d$TIME, -.d$EVID), ]
}

.edT5Model <- function() {
  ini({
    # DISPLACED from the simulated truth (1.63, 1.55, -2.4, -2.4, rho 0.5) on
    # purpose: starting a parameter at its true value makes "the fit recovers
    # it" pass without the estimator moving anything, which is the same trap
    # bWT falls into below.
    lclm <- 1.9; lv1m <- 1.8
    lclrv <- -2.0; lv1rv <- -2.0
    bWT <- 0
    eta.cl + eta.v1 ~ c(1, 0.3, 1)
    dist(eta.cl) ~ dgamma(shape = 1 / exp(lclrv),
                          rate = 1 / (exp(lclrv) *
                                      exp(lclm + bWT * log(WT / 70))))
    dist(eta.v1) ~ dgamma(shape = 1 / exp(lv1rv),
                          rate = 1 / (exp(lv1rv) * exp(lv1m)))
    prop.sd <- 0.1
  })
  model({
    cl <- eta.cl; v <- eta.v1
    linCmt() ~ prop(prop.sd)
  })
}

# A variant with the residual variance bounded from below.  `prop()` alone makes
# the objective DISCONTINUOUS in the coefficient on this arm: a subject's inner
# MAP switches mode, its prediction collapses to ~6e-6, and since the
# proportional variance is (prop.sd*IPRED)^2 the log-variance term REWARDS the
# collapse -- objf 20370 at bWT 0.09 against 25278 at 0.10, on identical
# records.  add.sd is FIXED rather than estimated because a free additive term
# is driven to 0 and the degeneracy comes straight back.
.edT5ModelBounded <- function() {
  ini({
    lclm <- 1.9; lv1m <- 1.8
    lclrv <- -2.0; lv1rv <- -2.0
    bWT <- 0.1
    eta.cl + eta.v1 ~ c(1, 0.3, 1)
    dist(eta.cl) ~ dgamma(shape = 1 / exp(lclrv),
                          rate = 1 / (exp(lclrv) *
                                      exp(lclm + bWT * log(WT / 70))))
    dist(eta.v1) ~ dgamma(shape = 1 / exp(lv1rv),
                          rate = 1 / (exp(lv1rv) * exp(lv1m)))
    add.sd <- fix(0.01); prop.sd <- 0.1
  })
  model({
    cl <- eta.cl; v <- eta.v1
    linCmt() ~ add(add.sd) + prop(prop.sd)
  })
}

nmTest({
  test_that("T5: the degenerate arm fits, and prop.sd stays inside its bound", {
    .d <- .edT5Data()
    .f <- suppressMessages(suppressWarnings(
      nlmixr2(.edT5Model(), .d, est = "focei",
              control = foceiControl(print = 0L, covMethod = ""))))
    .p <- setNames(.f$parFixedDf$Estimate, rownames(.f$parFixedDf))

    # prop.sd carries lower = 0 and a true value of 0.10.  An earlier version of
    # this design drove it to -0.889 and then failed on a bounds violation; a
    # bounded parameter leaving its bound is a defect whatever else happens, so
    # this is asserted on its own.
    expect_true(.p[["prop.sd"]] > 0)
    expect_equal(.p[["prop.sd"]], 0.10, tolerance = 0.5)

    # The structural parameters are started away from truth and must come back
    # to it, so these cannot pass on a frozen fit.  Measured at 120 subjects:
    # lclm 1.9 -> 1.6228, lv1m 1.8 -> 1.5440, lclrv -2.0 -> -2.2394,
    # rxCor 0.3 -> 0.6017, against truth 1.63, 1.55, -2.4 and 0.5.
    expect_equal(.p[["lclm"]], 1.63, tolerance = 0.15)
    expect_equal(.p[["lv1m"]], 1.55, tolerance = 0.15)
    # directional rather than tight: the relative variances move toward truth
    # from -2.0 but do not arrive, and the copula moves well off its 0.3 start
    expect_true(.p[["lclrv"]] < -2.05)
    expect_true(.p[["rxCor.eta.v1.eta.cl"]] > 0.4)

    # NOTE: bWT is deliberately NOT asserted here.  Its start and its truth are
    # both 0, so an estimator that never moves it "recovers" it perfectly --
    # which is exactly what happens today (see the skipped test below).  Reading
    # bWT ~ 0 on this arm as evidence that the covariate machinery works is the
    # trap this file exists to avoid.
  })

  test_that("a covariate coefficient IS recovered, given a non-zero slope start", {
    # The arm with a real effect.  Two things have to be right for this to work,
    # and both were mistaken for the feature being broken:
    #
    #   1. the slope must not start at exactly 0.  Measured on this arm, focei
    #      returns 0.0017 from a start of 0 and 0.7353 from a start of 0.1 --
    #      the objective's central difference at 0 is -53, so the gradient is
    #      there; from 0 the search takes one ~0.001 step, the objective change
    #      falls under tolerance, and it declares convergence.  A coefficient
    #      conventionally starts at 0, which is the one value that traps it.
    #   2. the residual variance must be bounded below (see
    #      .edT5ModelBounded), or the objective is not even continuous in the
    #      coefficient.
    .d <- .edT5Data(bWT = 0.75)
    .f <- suppressMessages(suppressWarnings(
      nlmixr2(.edT5ModelBounded(), .d, est = "focei",
              control = foceiControl(print = 0L, covMethod = ""))))
    .p <- setNames(.f$parFixedDf$Estimate, rownames(.f$parFixedDf))
    # recovery, not merely movement
    expect_equal(.p[["bWT"]], 0.75, tolerance = 0.25)
    # and it must have left the trapped region entirely
    expect_true(.p[["bWT"]] > 0.3)
  })

  test_that("T4: an observation-count imbalance does not move the coefficient", {
    # The record-weighting assertion.  The declared M-step weights a record by
    # 1/n_i so each SUBJECT contributes one unit however often it was observed;
    # a per-record weighting would instead let the observation count pull the
    # coefficient.  Choosing the weighting proves nothing -- this measures it,
    # by fitting the SAME subjects and the same etas twice, once as simulated
    # and once thinned so the obs/subject ratio roughly doubles.
    #
    # Measured on the full 120-subject arm: balanced bWT 0.4088 (SE 0.1342) at
    # 1.8:1, thinned 0.3784 (SE 0.4464) at 3.5:1 -- a difference of 0.0305, or
    # 0.07 pooled standard errors.
    .d <- .edT5Data(bWT = 0.75)
    # thin by POSITION within subject, after the assay limit, so the retained
    # records are not selected by time (which would confound the comparison)
    .keep <- unlist(lapply(split(seq_len(nrow(.d)), .d$ID), function(.i) {
      .obs <- .i[.d$EVID[.i] == 0]
      .dose <- .i[.d$EVID[.i] != 0]
      if (length(.obs) <= 2L) return(.i)
      c(.dose, .obs[seq(1L, length(.obs), by = 2L)])
    }), use.names = FALSE)
    .thin <- .d[sort(.keep), ]

    .fit <- function(.dat) {
      .f <- suppressMessages(suppressWarnings(
        nlmixr2(.edT5ModelBounded(), .dat, est = "focei",
                control = foceiControl(print = 0L, covMethod = ""))))
      unname(.f$parFixedDf["bWT", "Estimate"])
    }
    .b1 <- .fit(.d)
    .b2 <- .fit(.thin)

    # both must have escaped the trapped region, or "they agree" is vacuous --
    # two fits both stuck near 0 would also "agree"
    expect_true(.b1 > 0.2)
    expect_true(.b2 > 0.2)
    # and then they must actually agree despite the different record counts
    expect_equal(.b2, .b1, tolerance = 0.35)
  })

  test_that("the coefficient is at least IN the outer problem, not frozen", {
    # Weaker than recovery but still worth pinning: a coefficient that never
    # moves at all is a different (and worse) defect than one whose search
    # stalls, and the two are indistinguishable on the degenerate arm.
    .d <- .edT5Data()
    .m <- .edT5Model()
    .m <- rxode2::ini(.m, bWT = 0.5)
    .f <- suppressMessages(suppressWarnings(
      nlmixr2(.m, .d, est = "focei",
              control = foceiControl(print = 0L, covMethod = ""))))
    .b <- .f$parFixedDf$Estimate[match("bWT", rownames(.f$parFixedDf))]
    expect_false(isTRUE(all.equal(.b, 0.5, tolerance = 1e-6)))
  })
})

# .etaDistAddCovariate(): adding a covariate to one ROLE of one declaration
# (plan phase 3.4).  rxode2's expansion hoists each family argument onto its own
# rxEdA.<eta>.<role> line, so this only has to add a term at the source -- the
# declaration -- and the anchor picks it up.
nmTest({
  .edcBase <- function() {
    ini({
      lclm <- 1.63; lclrv <- -2.4; prop.sd <- 0.1
      eta.cl ~ 1
      dist(eta.cl) ~ dgamma(shape = 1 / exp(lclrv),
                            rate = 1 / (exp(lclrv) * exp(lclm)))
    })
    model({
      cl <- eta.cl; v <- 5
      linCmt() ~ prop(prop.sd)
    })
  }

  test_that("a covariate is added to the requested role's anchor", {
    .u <- nlmixr2est:::.etaDistAddCovariate(
      suppressMessages(nlmixr2est::nlmixr2(.edcBase)),
      "eta.cl", "rate", "WT", shape = "power", center = 70)
    .ln <- vapply(.u$lstExpr, function(.z) paste(deparse(.z), collapse = " "),
                  character(1))
    .rate <- .ln[grepl("^rxEdA[.]eta[.]cl[.]rate", .ln)]
    .shape <- .ln[grepl("^rxEdA[.]eta[.]cl[.]shape", .ln)]
    # exactly one of each: re-expanding a ui that already carries the expansion
    # would give two sets of anchors and two decoders for one eta
    expect_length(.rate, 1L)
    expect_length(.shape, 1L)
    expect_length(.ln[grepl("^eta.cl <- gammapInv", .ln)], 1L)
    # the covariate is on the RATE and nowhere else
    expect_true(grepl("WT", .rate, fixed = TRUE))
    expect_false(grepl("WT", .shape, fixed = TRUE))
    # the coefficient exists and does NOT start at 0 -- a slope started at
    # exactly 0 has no magnitude for the outer search to scale by
    .b <- .u$iniDf[.u$iniDf$name == "beta.eta.cl.rate.WT", ]
    expect_equal(nrow(.b), 1L)
    expect_true(.b$est != 0)
  })

  test_that(".etaDistAddCovariate refuses what it cannot do", {
    .u <- suppressMessages(nlmixr2est::nlmixr2(.edcBase))
    # a role the family does not have
    expect_error(nlmixr2est:::.etaDistAddCovariate(.u, "eta.cl", "df", "WT"),
                 "has no role")
    # not a declared random effect
    expect_error(nlmixr2est:::.etaDistAddCovariate(.u, "eta.nope", "rate", "WT"),
                 "not a declared random effect")
    # the same covariate twice on the same role
    .u2 <- nlmixr2est:::.etaDistAddCovariate(.u, "eta.cl", "rate", "WT",
                                             center = 70)
    expect_error(nlmixr2est:::.etaDistAddCovariate(.u2, "eta.cl", "rate", "WT",
                                                   center = 70),
                 "already")
  })

  test_that("an added covariate is estimated, with the role's sign", {
    # The arm's true effect is +0.75 on the gamma MEAN.  This adds the term to
    # the RATE, and rate = 1/(rv*mean), so the same effect reads as -0.75 there
    # -- which is why `rate` is a separate role from `scale` rather than folded
    # into it.  Measured on the full 120-subject arm: -0.6010 against the
    # hand-written mean-scale model's +0.6020, same objective (157.79).
    .d <- .edT5Data(bWT = 0.75)
    .u <- nlmixr2est:::.etaDistAddCovariate(
      suppressMessages(nlmixr2est::nlmixr2(.edT5ModelBounded())),
      "eta.v1", "rate", "WT", shape = "power", center = 70)
    .f <- suppressMessages(suppressWarnings(
      nlmixr2(.u, .d, est = "focei",
              control = foceiControl(print = 0L, covMethod = ""))))
    expect_true("beta.eta.v1.rate.WT" %in% rownames(.f$parFixedDf))
    # it moved off its 0.1 start rather than sitting there
    .b <- unname(.f$parFixedDf["beta.eta.v1.rate.WT", "Estimate"])
    expect_false(isTRUE(all.equal(.b, 0.1, tolerance = 1e-6)))
  })
})
