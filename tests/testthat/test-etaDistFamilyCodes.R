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

})
