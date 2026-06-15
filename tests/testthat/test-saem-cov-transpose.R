nmTest({
  test_that("SAEM linFim covariance with a single estimated covariate parameter (drop=FALSE regression)", {
    # Regression for the wrong-transpose bug in calc.COV() (R/saem_fit_aux.R).
    #
    # The covariance-via-linearization code multiplies DFi %*% t(Ai[cov.est.ix, ]).
    # `cov.est.ix` selects the estimated covariate-model parameters; the intended
    # orientation of t(Ai[cov.est.ix, ]) is nphi x npar so that DFi (nobs x nphi)
    # %*% (.) gives one design column per estimated parameter (nobs x npar).
    #
    # When exactly ONE parameter is estimated, Ai[cov.est.ix, ] dropped to a
    # length-nphi vector, so t() produced a 1 x nphi row instead of nphi x 1.  For
    # nphi > 1 that is non-conformable: calc.COV() errored, and end-to-end the SAEM
    # fit failed in the fallback path ("Not a matrix") instead of returning the
    # linearized-FIM covariance.  drop = FALSE keeps it nphi x npar in all cases.

    d <- theo_sd
    d$logWt <- log(d$WT / 70)
    d <- d[d$ID %in% 1:6, ]

    # nphi = 3 (ka, cl, v); both typical values fixed so the ONLY estimated entry
    # of the covariate-parameter matrix is the WT-on-CL slope -> sum(cov.est.ix)==1.
    one.cmt.cov <- function() {
      ini({
        tka <- fix(log(1.5))
        tcl <- fix(log(2.7))
        tv  <- fix(log(31))
        cl.wt <- 0.75
        eta.ka ~ 0.6
        eta.cl ~ 0.3
        eta.v  ~ 0.1
        add.sd <- 0.7
      })
      model({
        ka <- exp(tka + eta.ka)
        cl <- exp(tcl + cl.wt * logWt + eta.cl)
        v  <- exp(tv + eta.v)
        d/dt(depot)  <- -ka * depot
        d/dt(center) <-  ka * depot - (cl / v) * center
        cp <- center / v
        cp ~ add(add.sd)
      })
    }

    ctl <- saemControl(nBurn = 30, nEm = 30, nmc = 3, print = 0, seed = 1L,
                       covMethod = "linFim")

    # Before the fix this throws inside the covariance step; after it completes.
    fit <- .nlmixr(one.cmt.cov, d, est = "saem", control = ctl)

    # End-to-end: the linearized FIM covariance was actually produced (not a
    # fallback), and it is the expected 1 x 1 covariance of the single slope.
    expect_equal(fit$covMethod, "linFim")
    expect_equal(dim(fit$cov), c(1L, 1L))
    expect_true(is.finite(fit$cov[[1L]]))
    expect_gt(fit$cov[[1L]], 0)

    # Value/orientation check.  The (k, k) entry of the linearized FIM (= X'X) is
    # the squared norm of that parameter's design column and is invariant to which
    # other parameters are estimated.  So the single-parameter FIM (1 / variance)
    # must equal the cl.wt diagonal of the FIM computed when every parameter is
    # treated as estimated -- this pins the correct value, not just "no error".
    .saem <- fit$env$saem
    .cfg  <- attr(.saem, "saem.cfg")
    .nphi <- .cfg$nphi1 + .cfg$nphi0
    .covEstIx <- function(cfg) {
      .fixed <- c(matrix(names(cfg$inits$theta), ncol = .nphi, byrow = TRUE)) == "FIXED"
      .present <- !is.na(c(matrix(cfg$inits$theta, ncol = .nphi, byrow = TRUE)))
      !.fixed & .present
    }

    .ixSingle <- .covEstIx(.cfg)
    expect_equal(sum(.ixSingle), 1L)              # the reproduction condition

    .covSingle <- suppressMessages(nlmixr2est:::calc.COV(.saem))  # errors on buggy code
    expect_equal(dim(.covSingle), c(1L, 1L))

    # Same fitted object, but pretend the (fixed) typical values are estimated too.
    .cfgAll <- .cfg
    names(.cfgAll$inits$theta)[names(.cfgAll$inits$theta) == "FIXED"] <- ""
    .saemAll <- .saem
    attr(.saemAll, "saem.cfg") <- .cfgAll
    .ixAll <- .covEstIx(.cfgAll)
    expect_gt(sum(.ixAll), 1L)                     # multi-param path (no drop)

    .idx <- match(which(.ixSingle), which(.ixAll)) # cl.wt column in the multi cov
    .covAll <- suppressMessages(nlmixr2est:::calc.COV(.saemAll))
    .fimAll <- solve(.covAll)

    expect_equal(unname(1 / .covSingle[[1L]]),
                 unname(.fimAll[.idx, .idx]),
                 tolerance = 1e-6)
  })
})
