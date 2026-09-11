nmTest({
  # fit$cor assumed $cov was a matrix, so a fit without a covariance
  # (covMethod="") errored instead of returning NULL, which killed print()
  # on a fit reloaded from a saved copy (#1038).
  one.cmt <- function() {
    ini({
      tka <- 0.45
      tcl <- 1
      tv <- 3.45
      eta.cl ~ 0.09
      add.sd <- 0.7
    })
    model({
      ka <- exp(tka)
      cl <- exp(tcl + eta.cl)
      v <- exp(tv)
      d/dt(depot) <- -ka * depot
      d/dt(center) <- ka * depot - cl / v * center
      cp <- center / v
      cp ~ add(add.sd)
    })
  }

  # nlmixr2save reloads a fit with list2env(), which parents the fit
  # environment on globalenv() instead of emptyenv(); there "cov" and "cor"
  # resolve to stats::cov/stats::cor unless the lookup is local.
  .asReloaded <- function(fit) {
    .env <- fit$env
    .re <- list2env(mget(ls(.env, all.names = TRUE), envir = .env),
      envir = new.env(parent = globalenv())
    )
    class(.re) <- c("nlmixr2FitCore", paste0("nlmixr2.", fit$est))
    .re
  }

  test_that("$cor is NULL (not an error) when covMethod='' (#1038)", {
    fit <- .nlmixr(one.cmt, nlmixr2data::theo_sd,
      est = "focei",
      control = foceiControl(
        print = 0, maxInnerIterations = 1, maxOuterIterations = 1,
        eval.max = 1, covMethod = ""
      )
    )

    expect_null(fit$cov)
    expect_null(fit$cor)
    expect_error(capture.output(print(fit)), NA)

    .re <- .asReloaded(fit)
    # without the local lookup these pick up stats::cov / stats::cor
    expect_null(.re$cov)
    expect_null(.re$cor)
    expect_error(capture.output(print(.re)), NA)
    expect_false(any(grepl("\\$cor", capture.output(print(.re)))))
  })

  test_that("$cor is still the theta correlation when a covariance exists (#1038)", {
    fit <- .nlmixr(one.cmt, nlmixr2data::theo_sd,
      est = "focei",
      control = foceiControl(
        print = 0, maxInnerIterations = 1, maxOuterIterations = 1,
        eval.max = 1, covMethod = "r"
      )
    )
    skip_if(is.null(fit$cov), "no covariance was calculated")

    .cor <- fit$cor
    expect_true(is.matrix(.cor))
    expect_equal(dimnames(.cor), dimnames(fit$cov))
    # diagonal is the standard error, off-diagonal the correlation
    expect_equal(diag(.cor), sqrt(diag(fit$cov)))
    expect_equal(.cor[lower.tri(.cor)],
      stats::cov2cor(fit$cov)[lower.tri(.cor)]
    )

    .re <- .asReloaded(fit)
    expect_equal(.re$cor, .cor)
    # the correlation line is reachable again; it was gated on
    # exists("cor", x$env), which is never true for a locally fit model
    expect_true(any(grepl("\\$cor", capture.output(print(.re)))))
  })
})
