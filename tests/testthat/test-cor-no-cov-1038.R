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
  # `drop` mimics the items nlmixr2save leaves out of the saved fit ("model")
  .asReloaded <- function(fit, parent = globalenv(), drop = character(0)) {
    .env <- fit$env
    .nm <- setdiff(ls(.env, all.names = TRUE), drop)
    .re <- list2env(mget(.nm, envir = .env),
      envir = new.env(parent = parent)
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
    expect_false(any(grepl("\\$cor", capture.output(print(fit)))))

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

    # The correlation line is reachable again.  It has to be asserted on the
    # LOCAL fit: the gate it replaced, exists("cor", x$env), is never true for
    # a locally fit model but IS true for a reloaded one (it finds stats::cor),
    # so a reloaded-only assertion passes with the fix reverted.
    expect_true(any(grepl("\\$cor", capture.output(print(fit)))))

    .re <- .asReloaded(fit)
    expect_equal(.re$cor, .cor)
    expect_true(any(grepl("\\$cor", capture.output(print(.re)))))

    # the strong-correlation branch calls .getCorPrint(), which was unreachable
    # for a local fit before this change
    withr::local_options(list(nlmixr2.strong.corr = 0))
    .out <- capture.output(print(fit))
    expect_true(any(grepl("strong fixed parameter correlations", .out)))
  })

  test_that("$cor keeps a zero-variance row out of cov2cor (#1038)", {
    .env <- new.env(parent = emptyenv())
    .nm <- list(c("a", "b", "c"), c("a", "b", "c"))
    assign("cov", matrix(c(
      4, 1, 0,
      1, 9, 0,
      0, 0, 0
    ), 3, 3, dimnames = .nm), envir = .env)
    .lst <- list(.env, FALSE)
    class(.lst) <- c("cor", "nmObjGet")

    expect_warning(.cor <- nmObjGet(.lst), NA)
    expect_equal(diag(.cor), c(a = 2, b = 3, c = 0))
    expect_equal(.cor["a", "b"], 1 / 6)
    expect_true(all(is.na(.cor["c", c("a", "b")])))
    # .getR() drops the NA row, so print() sees only the real correlation
    expect_equal(unname(.getR(.cor)), 1 / 6)
  })

  test_that("$cor is NULL for a non-matrix $cov (#1038)", {
    .lst <- list(new.env(parent = emptyenv()), FALSE)
    class(.lst) <- c("cor", "nmObjGet")
    for (.v in list(NULL, stats::cov, "a", data.frame(a = 1), matrix(1:6, 2, 3))) {
      assign("cov", .v, envir = .lst[[1]])
      expect_null(nmObjGet(.lst))
    }
    # an empty square matrix is a covariance, so $cor mirrors it like $cov does
    assign("cov", matrix(numeric(0), 0, 0), envir = .lst[[1]])
    expect_equal(nmObjGet(.lst), matrix(numeric(0), 0, 0))
  })

  test_that("a reloaded fit does not read fit items out of its parent (#1038)", {
    fit <- .nlmixr(one.cmt, nlmixr2data::theo_sd,
      est = "focei",
      control = foceiControl(
        print = 0, maxInnerIterations = 1, maxOuterIterations = 1,
        eval.max = 1, covMethod = "r"
      )
    )

    # stands in for the user workspace a reloaded fit is parented on
    .shadow <- new.env(parent = emptyenv())
    assign("cov", matrix(1, 1, 1, dimnames = list("bogus", "bogus")), envir = .shadow)
    assign("covList", list(bogus = matrix(1, 1, 1)), envir = .shadow)
    assign("ranef", "bogus", envir = .shadow)
    assign("mixNum", "bogus", envir = .shadow)
    assign("parHistData", "bogus", envir = .shadow)

    .re <- .asReloaded(fit, parent = .shadow,
      drop = c("cov", "ranef", "mixNum", "parHistData")
    )
    expect_null(.re$cov)
    expect_null(.re$cor)
    expect_null(.re$ranef)
    expect_null(.re$mixNum)
    expect_null(.re$parHist)
    expect_false(any(grepl("bogus", capture.output(print(.re)))))
  })
})
