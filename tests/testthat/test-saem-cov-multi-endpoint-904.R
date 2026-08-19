# calc.COV's covFull residual variance block used to build dVi/d(residual param)
# from every observation row of the subject, across ALL endpoints, instead of only
# the endpoint that residual parameter belongs to (#904).  Two unrelated endpoints
# (sharing no theta/eta) make this directly checkable: masked to its own endpoint,
# a residual parameter's SE from the joint fit must reproduce the SE from fitting
# that endpoint alone, since nothing about the other endpoint's data or parameters
# should be able to leak in.  Multi-iteration fits -- weekly batch.

nmTest({
  test_that("saem covFull residual variance is masked to its own endpoint (#904)", {
    skip_on_cran()
    # pin threads: the simulated data and the fits themselves must be
    # reproducible regardless of the runner's core count
    .oldThreads <- rxode2::getRxThreads()
    on.exit(rxode2::setRxThreads(.oldThreads), add = TRUE)
    rxode2::setRxThreads(1L)

    twoEp <- function() {
      ini({
        tka <- 0.45; tcl <- 1; tv <- 3.45
        eta.ka ~ 0.6; eta.cl ~ 0.3; eta.v ~ 0.1
        add.sd <- 0.7
        tbio <- 5
        eta.bio ~ 0.25
        pdadd.sd <- 0.15
      })
      model({
        ka <- exp(tka + eta.ka); cl <- exp(tcl + eta.cl); v <- exp(tv + eta.v)
        d/dt(depot) <- -ka * depot
        d/dt(center) <- ka * depot - cl / v * center
        cp <- center / v
        bio <- tbio + eta.bio                # unrelated endpoint: no shared theta/eta with cp
        cp ~ add(add.sd)
        # deliberately a DIFFERENT residual type (proportional, not additive) than cp's:
        # with matching types the unmasked pre-fix dVi/da columns for both parameters
        # were numerically identical, making the old code's blocB exactly singular and
        # solve() fail -- so the fix only showed up as a missing row, not a wrong value.
        # A mismatched type makes the pre-fix result wrong-but-present, so this test
        # actually exercises the SE comparison below rather than an earlier NULL guard.
        bio ~ prop(pdadd.sd) | bio
      })
    }
    cpOnly <- function() {
      ini({
        tka <- 0.45; tcl <- 1; tv <- 3.45
        eta.ka ~ 0.6; eta.cl ~ 0.3; eta.v ~ 0.1
        add.sd <- 0.7
      })
      model({
        ka <- exp(tka + eta.ka); cl <- exp(tcl + eta.cl); v <- exp(tv + eta.v)
        d/dt(depot) <- -ka * depot
        d/dt(center) <- ka * depot - cl / v * center
        cp <- center / v
        cp ~ add(add.sd)
      })
    }
    bioOnly <- function() {
      ini({ tbio <- 5; eta.bio ~ 0.25; pdadd.sd <- 0.15 })
      model({ bio <- tbio + eta.bio; bio ~ prop(pdadd.sd) })
    }

    .testSeed(2); rxode2::rxSetSeed(2)
    .N <- 30
    .ev <- rxode2::et(amt = 100, cmt = "depot", id = 1:.N)
    .ev <- rxode2::et(.ev, seq(0.5, 24, by = 3), cmt = "cp")
    .ev <- rxode2::et(.ev, seq(0.5, 24, by = 3), cmt = "bio")
    .d <- as.data.frame(rxode2::rxSolve(twoEp, .ev, addDosing = TRUE))
    .dose <- .d[.d$evid != 0, c("id", "time", "CMT", "amt", "evid")]
    .dose$dv <- NA_real_
    .obs <- .d[.d$evid == 0, c("id", "time", "CMT", "sim")]
    .obs$amt <- 0; .obs$evid <- 0
    names(.obs)[names(.obs) == "sim"] <- "dv"
    .dat <- rbind(.dose[, c("id", "time", "dv", "CMT", "amt", "evid")],
                  .obs[, c("id", "time", "dv", "CMT", "amt", "evid")])
    .dat <- .dat[order(.dat$id, .dat$time, -.dat$evid), ]

    # dosing + cp observations are exactly the rows the cp-only model sees -- "the
    # first endpoint's data is unchanged" (issue #904's verification recipe)
    .datCp <- .dat[.dat$CMT %in% c(1L, 3L), ]
    .datBio <- .dat[.dat$CMT == 4L, ]

    .ctl <- saemControl(nBurn = 150, nEm = 200, print = 0, seed = 1L,
                        calcTables = FALSE, covMethod = "linFim")
    .fJoint <- .nlmixr(twoEp, .dat, est = "saem", control = .ctl)
    .fCp    <- .nlmixr(cpOnly, .datCp, est = "saem", control = .ctl)
    .fBio   <- .nlmixr(bioOnly, .datBio, est = "saem", control = .ctl)

    .getVarCov <- function(f) {
      .s <- f$saem
      attr(.s, "env") <- f$env
      attr(suppressWarnings(calc.COV(.s)), "varCov")
    }
    .vcJoint <- .getVarCov(.fJoint)
    .vcCp    <- .getVarCov(.fCp)
    .vcBio   <- .getVarCov(.fBio)

    expect_true(all(c("add.sd", "pdadd.sd") %in% rownames(.vcJoint)))
    expect_true(all(is.finite(diag(.vcJoint))))
    expect_true(all(diag(.vcJoint) > 0))

    # the mechanism: masked to its own endpoint, each residual parameter's SE from
    # the joint fit reproduces the single-endpoint reference (tight tolerance --
    # nothing about the OTHER endpoint should be able to move this at all).  These
    # SEs are ~0.01 in magnitude, well under a naive absolute tolerance -- compare
    # the RATIO to 1 so the check is genuinely relative (expect_equal(x, y, tolerance=)
    # falls back to an absolute difference once both values are small, which would
    # make a tolerance chosen for O(1) numbers pass almost regardless of x vs y).
    expect_equal(sqrt(.vcJoint["add.sd", "add.sd"]) / sqrt(.vcCp["add.sd", "add.sd"]),
                 1, tolerance = 0.15)
    expect_equal(sqrt(.vcJoint["pdadd.sd", "pdadd.sd"]) / sqrt(.vcBio["pdadd.sd", "pdadd.sd"]),
                 1, tolerance = 0.15)
  })
})
