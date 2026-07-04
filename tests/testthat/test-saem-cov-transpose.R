nmTest({
  test_that("SAEM linFim covariance with a single estimated covariate parameter (drop=FALSE regression)", {
    # Regression for the wrong-transpose bug in calc.COV() (R/saem_fit_aux.R)
    # when exactly one covariate-model parameter is estimated.

    d <- theo_sd
    d$logWt <- log(d$WT / 70)
    d <- d[d$ID %in% 1:6, ]

    # nphi = 3 (ka, cl, v); tka/tv fixed so cl.wt is the only estimated covariate parameter
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

    fit <- .nlmixr(one.cmt.cov, d, est = "saem", control = ctl)

    # linFim covariance produced end-to-end (not a fallback)
    expect_equal(fit$covMethod, "linFim")
    expect_equal(dim(fit$cov), c(1L, 1L))
    expect_true(is.finite(fit$cov[[1L]]))
    expect_gt(fit$cov[[1L]], 0)

    # value check: single-parameter FIM (1/variance) must equal the cl.wt
    # diagonal of the FIM when every parameter is treated as estimated
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
