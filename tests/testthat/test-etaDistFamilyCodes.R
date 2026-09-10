# `.etaDistFamilyCode()` returns the ROW INDEX of the family in
# `lotri::lotriEtaDists()`, and that integer is what the C++ dispatch switches
# on -- `src/etaDistFam.h`'s `#define RXETADIST_GAMMA 13` means "the 13th row of
# lotri's table", nothing more.
#
# So inserting or reordering a family in lotri's `.lotriEtaDistDefs` silently
# rebinds EVERY family at or after the insertion point: a model declaring
# dgamma() would be fitted as whatever now sits at row 13, with no error from
# either package.  Adding a COLUMN to that table is safe; adding a ROW is not.
#
# This pins the whole binding so such a change breaks a test instead of a fit.
# If a family genuinely needs to be added, append it at the END, extend this
# table, and add the matching `#define` -- in that order.

nmTest({

  test_that("the family name -> code binding matches the C++ catalog", {
    skip_if_not_installed("lotri")
    .expected <- c(
      dnorm              =  1L,   # RXETADIST_NORM        mean, sd
      stdNormal          =  2L,   # RXETADIST_STDNORMAL   (none)
      studentT           =  3L,   # RXETADIST_STUDENTT    nu, mu, sigma
      dcauchy            =  4L,   # RXETADIST_CAUCHY      location, scale
      doubleExponential  =  5L,   # RXETADIST_DBLEXP      mu, sigma
      dlogis             =  6L,   # RXETADIST_LOGIS       location, scale
      gumbel             =  7L,   # RXETADIST_GUMBEL      mu, beta
      dlnorm             =  8L,   # RXETADIST_LNORM       meanlog, sdlog
      dchisq             =  9L,   # RXETADIST_CHISQ       df
      invChiSquare       = 10L,   # RXETADIST_INVCHISQ    nu
      scaledInvChiSquare = 11L,   # RXETADIST_SCINVCHISQ  nu, sigma
      dexp               = 12L,   # RXETADIST_EXP         rate
      dgamma             = 13L,   # RXETADIST_GAMMA       shape, rate
      invGamma           = 14L,   # RXETADIST_INVGAMMA    alpha, beta
      dweibull           = 15L,   # RXETADIST_WEIBULL     shape, scale
      frechet            = 16L,   # RXETADIST_FRECHET     alpha, sigma
      rayleigh           = 17L,   # RXETADIST_RAYLEIGH    sigma
      pareto             = 18L,   # RXETADIST_PARETO      y_min, alpha
      paretoType2        = 19L,   # RXETADIST_PARETO2     mu, lambda, alpha
      dbeta              = 20L,   # RXETADIST_BETA        shape1, shape2
      betaProportion     = 21L,   # RXETADIST_BETAPROP    mu, kappa
      dunif              = 22L)   # RXETADIST_UNIF        min, max

    .tab <- lotri::lotriEtaDists()
    # a new family appended at the end is fine; a new one INSERTED is not
    expect_gte(nrow(.tab), length(.expected))
    expect_identical(.tab$name[seq_along(.expected)], names(.expected))

    for (.nm in names(.expected)) {
      expect_identical(nlmixr2est:::.etaDistFamilyCode(paste0(.nm, "()")),
                       .expected[[.nm]],
                       label = .nm)
    }
  })

  test_that("argument counts match what the C++ dispatch expects", {
    skip_if_not_installed("lotri")
    # rxEtaDistNarg() in src/etaDistFam.h is a second, hand-maintained copy of
    # this; if it and lotri disagree the M-step reads past the end of its
    # argument array or silently ignores an argument.
    .nReq <- c(dnorm = 2L, stdNormal = 0L, studentT = 3L, dcauchy = 2L,
               doubleExponential = 2L, dlogis = 2L, gumbel = 2L, dlnorm = 2L,
               dchisq = 1L, invChiSquare = 1L, scaledInvChiSquare = 2L,
               dexp = 1L, dgamma = 2L, invGamma = 2L, dweibull = 2L,
               frechet = 2L, rayleigh = 1L, pareto = 2L, paretoType2 = 3L,
               dbeta = 2L, betaProportion = 2L, dunif = 2L)
    .tab <- lotri::lotriEtaDists()
    for (.nm in names(.nReq)) {
      .w <- which(.tab$name == .nm)
      expect_length(.w, 1L)
      expect_identical(as.integer(.tab$nReq[.w]), .nReq[[.nm]], label = .nm)
    }
  })

  ## --------------------------------------------------------------- roles --

  test_that("roles come back per argument, in the family's own order", {
    skip_if_not_installed("lotri")
    skip_if(is.null(lotri::lotriEtaDists()$roles), "lotri has no roles column")
    expect_identical(nlmixr2est:::.etaDistRoles(quote(dgamma(shape = a, rate = b))),
                     c("shape", "rate"))
    ## the ORDER is the family's, not the call's -- the roles describe the
    ## catalog's positional argument list, which is what the M-step estimates
    expect_identical(nlmixr2est:::.etaDistRoles(quote(dgamma(rate = b, shape = a))),
                     c("shape", "rate"))
    expect_identical(nlmixr2est:::.etaDistRoles(quote(studentT(n, m, s))),
                     c("df", "location", "scale"))
    expect_identical(nlmixr2est:::.etaDistRoles(quote(dbeta(a, b))), c("shape1", "shape2"))
    expect_identical(nlmixr2est:::.etaDistRoles(quote(stdNormal())), character(0))
    ## an unknown family is "no role information", NOT "no roles"
    expect_identical(nlmixr2est:::.etaDistRoles(quote(dnotAFamily(a))), character(0))
    ## a string is accepted the same way the rest of this file's helpers do
    expect_identical(nlmixr2est:::.etaDistRoles("dexp(r)"), "rate")
  })

  test_that("roles line up one-for-one with the arguments the M-step fits", {
    skip_if_not_installed("lotri")
    .tab <- lotri::lotriEtaDists()
    skip_if(is.null(.tab$roles), "lotri has no roles column")
    ## The M-step indexes the declaration's arguments POSITIONALLY
    ## (`as.list(.cl)[-1]`), so a roles vector of a different length would
    ## attach a covariate to the wrong parameter rather than fail.
    for (.i in seq_len(nrow(.tab))) {
      .r <- strsplit(.tab$roles[.i], ",", fixed = TRUE)[[1]]
      .r <- .r[nzchar(.r)]
      expect_identical(length(.r), as.integer(.tab$nPar[.i]), label = .tab$name[.i])
    }
  })

  test_that("role groups are usable as covariate group keys", {
    skip_if_not_installed("lotri")
    skip_if(is.null(lotri::lotriEtaDists()$roles), "lotri has no roles column")
    .g <- nlmixr2est:::.etaDistRoleGroups(quote(dgamma(shape = a, rate = b)))
    expect_identical(names(.g), c("shape", "rate"))
    expect_identical(.g[["shape"]], 1L)
    expect_identical(.g[["rate"]], 2L)
    ## every argument lands in exactly one group -- a covariate attached to a
    ## role must reach one parameter, not zero and not two
    .n <- nlmixr2est:::.etaDistRoles(quote(paretoType2(m, l, a)))
    .g2 <- nlmixr2est:::.etaDistRoleGroups(quote(paretoType2(m, l, a)))
    expect_equal(sort(unlist(.g2, use.names = FALSE)), seq_along(.n))
    expect_identical(names(.g2), c("location", "scale", "shape"))
    ## unknown family -> NULL, so a caller declines rather than grouping
    ## every argument together under one empty key
    expect_null(nlmixr2est:::.etaDistRoleGroups(quote(dnotAFamily(a))))
  })

  test_that("the support endpoints are the roles that refuse covariates", {
    skip_if_not_installed("lotri")
    skip_if(is.null(lotri::lotriEtaDists()$roles), "lotri has no roles column")
    ## dunif's bounds and pareto's minimum ARE the support.  A subject-varying
    ## endpoint makes the density discontinuous in the parameter.
    expect_true(all(nlmixr2est:::.etaDistRoles(quote(dunif(lo, hi))) %in% nlmixr2est:::.etaDistRoleNoCovariate))
    expect_true("lower" %in% nlmixr2est:::.etaDistRoles(quote(pareto(ymin, alpha))))
    ## but pareto's shape is ordinary and must stay searchable
    expect_false("shape" %in% nlmixr2est:::.etaDistRoleNoCovariate)
    ## and nothing in the families the C++ M-step implements is a bound
    expect_false(any(nlmixr2est:::.etaDistRoles(quote(dgamma(a, b))) %in% nlmixr2est:::.etaDistRoleNoCovariate))
  })

})
