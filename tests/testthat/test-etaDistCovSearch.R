# Phase 3.4: the stepwise search over declared-distribution argument roles.
#
# The structural tests below are the ones that can fail cheaply and often.  The
# end-to-end search is guarded because it is many refits, but it is the only one
# that shows the driver SELECTS rather than merely runs -- the documented
# failure mode in this area is machinery that executes and decides nothing.

nmTest({
  test_that("a gamma declaration offers its shape and rate as targets", {
    .m <- function() {
      ini({
        lclm <- 1.63; lclrv <- -2.4
        eta.cl ~ 1
        dist(eta.cl) ~ dgamma(shape = 1 / exp(lclrv),
                              rate = 1 / (exp(lclrv) * exp(lclm)))
        prop.sd <- 0.1
      })
      model({ cl <- eta.cl; v <- 4.7; linCmt() ~ prop(prop.sd) })
    }
    .t <- .etaDistCovTargets(rxode2::rxode2(.m))
    expect_setequal(.t$role, c("shape", "rate"))
    expect_true(all(.t$eta == "eta.cl"))
    # the role is the group key, so each appears once
    expect_equal(anyDuplicated(.t$role), 0L)
  })

  test_that("a support endpoint is NOT offered as a target", {
    # a subject-varying bound makes the density discontinuous in the parameter,
    # so there is nothing for the M-step to follow -- this is the refusal, and
    # it must hold at the SEARCH boundary, not only inside the applier
    skip_if_not("lower" %in% .etaDistRoles(quote(dunif(a, b))))
    .m <- function() {
      ini({
        lo <- 0.5; hi <- 3
        eta.cl ~ 1
        dist(eta.cl) ~ dunif(lo, hi)
        prop.sd <- 0.1
      })
      model({ cl <- eta.cl; v <- 4.7; linCmt() ~ prop(prop.sd) })
    }
    .t <- .etaDistCovTargets(rxode2::rxode2(.m))
    expect_false(any(c("lower", "upper") %in% .t$role))
  })

  test_that("an undeclared model offers nothing, rather than everything", {
    .m <- function() {
      ini({ lcl <- 1.6; eta.cl ~ 0.1; prop.sd <- 0.1 })
      model({ cl <- exp(lcl + eta.cl); v <- 4.7; linCmt() ~ prop(prop.sd) })
    }
    expect_equal(NROW(.etaDistCovTargets(rxode2::rxode2(.m))), 0L)
    expect_error(etaDistCovarSearch(.m, data.frame(ID = 1)), "nothing to search")
  })
})

nmTest({
  .csModel <- function() {
    ini({
      lclm <- 1.63; lv1m <- 1.55
      lclrv <- -2.4; lv1rv <- -2.4
      eta.cl + eta.v1 ~ c(1, 0.5, 1)
      dist(eta.cl) ~ dgamma(shape = 1 / exp(lclrv),
                            rate = 1 / (exp(lclrv) * exp(lclm)))
      dist(eta.v1) ~ dgamma(shape = 1 / exp(lv1rv),
                            rate = 1 / (exp(lv1rv) * exp(lv1m)))
      prop.sd <- 0.1
    })
    model({ cl <- eta.cl; v <- eta.v1; linCmt() ~ prop(prop.sd) })
  }

  test_that("the candidate set is targets x covariates, one shape each", {
    .cov <- list(covNames = c("WT_p", "WT_l", "AGE_p"),
                 covRaw = c("WT", "WT", "AGE"),
                 covShape = c("power", "log", "power"),
                 covPop = c(70, 70, 40),
                 # covCanon is the existing "one column per covGroup" mask: the
                 # two WT shapes are one group, so only the first is canonical
                 covCanon = c(TRUE, FALSE, TRUE),
                 tvExcl = character(0))
    .c <- .etaDistCovCandidates(rxode2::rxode2(.csModel), .cov)
    # 2 etas x 2 roles x 2 canonical covariates
    expect_equal(NROW(.c), 8L)
    expect_setequal(unique(.c$cov), c("WT", "AGE"))
    # the non-canonical WT shape must NOT appear: two shapes of one covariate
    # are mutually exclusive, and proposing both lets them be selected together
    expect_false("log" %in% .c$shape)
  })

  test_that("a relation already in the model is not re-proposed", {
    .cov <- list(covNames = "WT", covRaw = "WT", covShape = "power",
                 covPop = 70, covCanon = TRUE, tvExcl = character(0))
    .all <- .etaDistCovCandidates(rxode2::rxode2(.csModel), .cov)
    expect_equal(NROW(.all), 4L)
    .less <- .etaDistCovCandidates(rxode2::rxode2(.csModel), .cov,
                                   already = .all[1, , drop = FALSE])
    expect_equal(NROW(.less), 3L)
    # and it is the RIGHT one that went
    expect_false(any(.less$eta == .all$eta[1] & .less$role == .all$role[1]))
  })

  test_that("an unscoreable candidate is reported, not silently ranked last", {
    # a relation that cannot be applied must travel with its reason: scoring it
    # as Inf would let it lose quietly and look considered
    .r <- .etaDistCovFit(rxode2::rxode2(.csModel), data.frame(ID = 1),
                         "focei", NULL,
                         data.frame(eta = "eta.cl", role = "rate",
                                    cov = "NOSUCHCOV", shape = "power",
                                    center = NA_real_,
                                    stringsAsFactors = FALSE))
    expect_false(.r$ok)
    expect_true(nzchar(.r$why))
  })
})

nmTest({
  test_that("the search picks the right covariate and eta out of noise", {
    skip_on_cran()
    # The only test here that shows the driver SELECTS.  A covariate search that
    # merely runs proves nothing: what has to be shown is that the true relation
    # beats a decoy on the same data, in one pass.
    #
    # Truth: WT on the cl declaration.  AGE is drawn independently and has no
    # effect at all, so it is the decoy; eta.v1 is the wrong-parameter decoy.
    #
    # Measured on this arm (60 subjects, focei, one forward step, ~41s):
    #
    #   eta.cl shape WT   dOFV -12.47   <- chosen
    #   eta.cl rate  WT   dOFV -11.61
    #   eta.v1 rate  WT   dOFV  -2.05
    #   eta.v1 shape WT   dOFV  -1.78
    #   eta.v1 rate  AGE  dOFV  -0.87
    #   eta.cl shape AGE  dOFV  -0.19
    #
    # so WT beats AGE by an order of magnitude and eta.cl beats eta.v1 by six
    # times, both far outside the noise.
    #
    # The ROLE is deliberately NOT asserted.  The truth is on the RATE and the
    # search picks SHAPE, by 0.85 objf units.  That is the data, not a defect:
    # a gamma's mean is shape/rate, so a covariate on either shifts it, and the
    # two differ only in the relative variance (rv = 1/shape), which rate does
    # not touch.  Asserting the role would pin a coin-flip and would also claim
    # a resolution this arm does not have -- see the plan's phase 3.6.
    .d <- .edT5Data(bWT = 0.75, timeVarying = FALSE)
    set.seed(7)
    .ids <- unique(.d$ID)
    .d <- merge(.d, data.frame(ID = .ids,
                               AGE = round(stats::rnorm(length(.ids), 40, 9))),
                by = "ID")
    .d <- .d[order(.d$ID, .d$TIME, -.d$EVID), ]
    .base <- function() {
      ini({
        lclm <- 1.63; lv1m <- 1.55
        lclrv <- -2.4; lv1rv <- -2.4
        eta.cl + eta.v1 ~ c(1, 0.5, 1)
        dist(eta.cl) ~ dgamma(shape = 1 / exp(lclrv),
                              rate = 1 / (exp(lclrv) * exp(lclm)))
        dist(eta.v1) ~ dgamma(shape = 1 / exp(lv1rv),
                              rate = 1 / (exp(lv1rv) * exp(lv1m)))
        prop.sd <- 0.1
      })
      model({ cl <- eta.cl; v <- eta.v1; linCmt() ~ prop(prop.sd) })
    }
    .r <- etaDistCovarSearch(.base, .d, est = "focei",
                             control = foceiControl(print = 0L, covMethod = ""),
                             maxSteps = 1L)
    # it selected SOMETHING -- a search that adds nothing would pass every
    # assertion phrased as "did not pick AGE"
    expect_equal(NROW(.r$relations), 1L)
    expect_identical(.r$relations$cov, "WT")
    expect_identical(.r$relations$eta, "eta.cl")
    # every candidate was scored, none silently dropped
    expect_equal(NROW(.r$path), 8L)
    expect_true(all(is.na(.r$path$why)))
    # and the margin is real, not a tie broken by ordering
    .wt <- min(.r$path$dOfv[.r$path$cov == "WT" & .r$path$eta == "eta.cl"])
    .age <- min(.r$path$dOfv[.r$path$cov == "AGE"])
    expect_true(.wt < .age - 5)
  })
})
