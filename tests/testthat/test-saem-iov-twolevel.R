# Unit tests for the two-level (Panhard & Samson) IOV handling in saem.  No
# fits here -- the fit-based tests live in test-saem-iov-twolevel-fit.R.

.twoLevelData <- function() {
  .d <- nlmixr2data::theo_md
  .d$occ <- 1
  .d$occ[.d$TIME >= 144] <- 2
  .d
}

# nolint start: object_usage_linter. rxode2 ini()/model() are NSE blocks:
# every assignment here is model specification, consumed by rxode2, not a
# local variable lintr can see used.
.twoLevelModel <- function() {
  ini({
    tka <- 0.45
    tcl <- 1
    tv <- 3.45
    add.sd <- 0.7
    eta.ka ~ 0.6
    eta.cl ~ 0.3
    eta.v ~ 0.1
    iov.cl ~ 0.1 | occ
  })
  model({
    ka <- exp(tka + eta.ka)
    cl <- exp(tcl + eta.cl + iov.cl)
    v <- exp(tv + eta.v)
    linCmt() ~ add(add.sd)
  })
}
# nolint end

test_that(".saemIovInfo describes the two-level structure", {
  .ui <- rxode2::rxUiDecompress(.twoLevelModel())
  .i <- .saemIovInfo(.ui, .twoLevelData())
  expect_true(is.list(.i))
  expect_equal(.i$occVar, "occ")
  expect_equal(.i$levels, c(1, 2))
  expect_equal(.i$pars$iov, "iov.cl")
  # the occasion term rides tcl, in the same additive position as eta.cl --
  # that is what makes it mu + b_i + c_ik
  expect_equal(.i$pars$theta, "tcl")
  expect_equal(.i$pars$eta, "eta.cl")
  expect_equal(.i$etaNames$iov.cl, c("rx.iov.cl.1", "rx.iov.cl.2"))

  # a model with no IOV at all
  # nolint start: object_usage_linter. rxode2 ini()/model() are NSE blocks:
  # every assignment here is model specification, consumed by rxode2, not a
  # local variable lintr can see used.
  .noIov <- function() {
    ini({
      tka <- 0.45
      tcl <- 1
      tv <- 3.45
      add.sd <- 0.7
      eta.ka ~ 0.6
    })
    model({
      ka <- exp(tka + eta.ka)
      cl <- exp(tcl)
      v <- exp(tv)
      linCmt() ~ add(add.sd)
    })
  }
  # nolint end
  expect_null(.saemIovInfo(rxode2::rxUiDecompress(.noIov()), .twoLevelData()))
})

test_that(".saemIovInfo declines what the two-level handling cannot take", {
  .d <- .twoLevelData()
  .d$occ2 <- ifelse(.d$TIME %% 2 < 1, 1, 2)
  # nolint start: object_usage_linter. rxode2 ini()/model() are NSE blocks:
  # every assignment here is model specification, consumed by rxode2, not a
  # local variable lintr can see used.
  .two <- function() {
    ini({
      tka <- 0.45
      tcl <- 1
      tv <- 3.45
      add.sd <- 0.7
      eta.ka ~ 0.6
      eta.cl ~ 0.3
      eta.v ~ 0.1
      iov.cl ~ 0.1 | occ
      iov.v ~ 0.1 | occ2
    })
    model({
      ka <- exp(tka + eta.ka)
      cl <- exp(tcl + eta.cl + iov.cl)
      v <- exp(tv + eta.v + iov.v)
      linCmt() ~ add(add.sd)
    })
  }
  # nolint end
  .why <- .saemIovInfo(rxode2::rxUiDecompress(.two()), .d)
  expect_true(is.character(.why))
  # the reason is collected into the fit's $runInfo, so it has to fit on a line
  expect_lt(nchar(.why), 75L)

  # a decline is what makes the shared rewrite run after all
  expect_equal(.iovNativeDecline(rxode2::rxUiDecompress(.two()), "saem", .d,
                                 saemControl(iovMethod = "twoLevel")),
               .why)
  # ... and no method without the attribute ever declines
  expect_null(.iovNativeDecline(rxode2::rxUiDecompress(.two()), "focei", .d,
                                foceiControl()))
})

test_that(".saemIovExpandUi writes the occasion term as a variance component", {
  .ui <- rxode2::rxUiDecompress(.twoLevelModel())
  .exp <- .saemIovExpandUi(.ui, .saemIovInfo(.ui, .twoLevelData()))
  .ini <- .exp$iniDf

  # no magnitude theta: the thetas are exactly the user's
  expect_equal(.ini$name[is.na(.ini$neta1)], c("tka", "tcl", "tv", "add.sd"))
  # one zero-mean eta per occasion level, carrying the declared variance
  expect_true(all(c("rx.iov.cl.1", "rx.iov.cl.2") %in% .ini$name))
  expect_equal(.ini$est[.ini$name == "rx.iov.cl.1"], 0.1)
  expect_equal(.ini$est[.ini$name == "rx.iov.cl.2"], 0.1)
  expect_true(all(.ini$condition[.ini$name %in% c("rx.iov.cl.1", "rx.iov.cl.2")] == "id"))
  expect_false("iov.cl" %in% .ini$name)

  # eta.cl stays mu-referenced to tcl -- the whole point of the additive form
  expect_true("eta.cl" %in% .exp$muRefDataFrame$eta)
  expect_equal(.exp$muRefDataFrame$theta[.exp$muRefDataFrame$eta == "eta.cl"], "tcl")
  # and the occasion etas get their own phi columns with means pinned at 0
  expect_equal(.exp$nonMuEtas, c("rx.iov.cl.1", "rx.iov.cl.2"))
  expect_true(all(.exp$saemFixed[c("rx.iov.cl.1", "rx.iov.cl.2")]))

  # the combination is on its own line: rxode2 refuses theta + eta1 + eta2
  expect_true(any(vapply(.exp$lstExpr,
                         function(l) identical(deparse1(l),
                                               "rx.iov.cl <- (occ == 1) * rx.iov.cl.1 + (occ == 2) * rx.iov.cl.2"),
                         logical(1))))
})

test_that("saemOmegaPool groups the occasion etas, and only those", {
  .ui <- rxode2::rxUiDecompress(.twoLevelModel())
  .exp <- .saemIovExpandUi(.ui, .saemIovInfo(.ui, .twoLevelData()))
  # saemEtaNames order is Gamma2_phi1 order, which is what the kernel indexes
  expect_equal(.exp$saemEtaNames,
               c("eta.ka", "eta.cl", "eta.v", "rx.iov.cl.1", "rx.iov.cl.2"))
  expect_equal(rxUiGet.saemOmegaPool(list(.exp)), c(0L, 0L, 0L, 1L, 1L))

  # the shared rewrite's line is a PRODUCT (magnitude * sum), so it must not be
  # picked up as a pool
  .legacy <- .uiApplyIov(.ui, "saem", .twoLevelData(),
                         saemControl(iovMethod = "theta"))$ui
  expect_true(all(rxUiGet.saemOmegaPool(list(.legacy)) == 0L))
})

test_that(".saemIovCollapseCov contracts the pooled occasion columns", {
  # the K per-occasion columns estimate ONE variance, so the covariance matrix
  # has to carry one row for them -- Var(mean(v_1..v_K)), not Var(v_1)
  .nm <- c("tka", "om.eta.ka", "om.rx.iov.cl.1", "om.rx.iov.cl.2")
  .cv <- matrix(c(4, 1, 2, 3,
                  1, 5, 6, 7,
                  2, 6, 10, -4,
                  3, 7, -4, 8), nrow = 4, byrow = TRUE,
                dimnames = list(.nm, .nm))
  .out <- .saemIovCollapseCov(.cv, list(iov.cl = c("rx.iov.cl.1", "rx.iov.cl.2")))

  expect_equal(rownames(.out), c("tka", "om.eta.ka", "om.iov.cl"))
  expect_equal(colnames(.out), c("tka", "om.eta.ka", "om.iov.cl"))
  # the untouched block is carried through unchanged
  expect_equal(.out["tka", "tka"], 4)
  expect_equal(.out["tka", "om.eta.ka"], 1)
  # Var(mean(v1,v2)) = (Var v1 + Var v2 + 2 Cov)/4 = (10 + 8 + 2*(-4))/4
  expect_equal(.out["om.iov.cl", "om.iov.cl"], (10 + 8 + 2 * -4) / 4)
  # Cov(mean(v1,v2), p) = (Cov(v1,p) + Cov(v2,p))/2
  expect_equal(.out["tka", "om.iov.cl"], (2 + 3) / 2)
  expect_equal(.out["om.eta.ka", "om.iov.cl"], (6 + 7) / 2)
  # symmetry survives
  expect_equal(.out, t(.out))

  # a matrix with nothing to pool, and an empty group list, are both untouched
  expect_identical(.saemIovCollapseCov(.cv, list()), .cv)
  expect_identical(.saemIovCollapseCov(.cv, list(iov.v = "rx.iov.v.1")), .cv)
})
