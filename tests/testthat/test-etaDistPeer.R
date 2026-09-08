## The M-step peer model's log densities (R/etaDistPeer.R).
##
## Two independent checks per family, neither of which goes through rxode2:
##
##   * against R's own log density, where one exists;
##   * exp(logDensity) integrates to 1 over the support.
##
## The second is what catches a wrong normalizing constant in the families
## written out by hand (gumbel, frechet, pareto, the inverse chi squares, ...),
## which have no R counterpart to compare against.  A constant that is wrong
## but does not depend on the parameters would still give the right MAXIMIZER,
## so a gradient check alone would not see it -- but it makes the reported
## objective meaningless and would break any comparison across families.

.edPeerSpec <- list(
  dnorm=list(p=c(mean=0.3, sd=1.4), lo=-Inf, hi=Inf,
             ref=function(x, p) stats::dnorm(x, p[1], p[2], log=TRUE)),
  studentT=list(p=c(nu=6, mu=0.2, sigma=1.3), lo=-Inf, hi=Inf,
                ref=function(x, p) stats::dt((x - p[2])/p[3], p[1], log=TRUE) - log(p[3])),
  dcauchy=list(p=c(location=0.1, scale=0.9), lo=-Inf, hi=Inf,
               ref=function(x, p) stats::dcauchy(x, p[1], p[2], log=TRUE)),
  doubleExponential=list(p=c(mu=0.2, sigma=1.1), lo=-Inf, hi=Inf, ref=NULL),
  dlogis=list(p=c(location=0.4, scale=1.2), lo=-Inf, hi=Inf,
              ref=function(x, p) stats::dlogis(x, p[1], p[2], log=TRUE)),
  gumbel=list(p=c(mu=0.3, beta=1.1), lo=-Inf, hi=Inf, ref=NULL),
  dlnorm=list(p=c(meanlog=0.2, sdlog=0.7), lo=0, hi=Inf,
              ref=function(x, p) stats::dlnorm(x, p[1], p[2], log=TRUE)),
  dchisq=list(p=c(df=4.5), lo=0, hi=Inf,
              ref=function(x, p) stats::dchisq(x, p[1], log=TRUE)),
  invChiSquare=list(p=c(nu=5), lo=0, hi=Inf, ref=NULL),
  scaledInvChiSquare=list(p=c(nu=5, sigma=1.3), lo=0, hi=Inf, ref=NULL),
  dexp=list(p=c(rate=1.7), lo=0, hi=Inf,
            ref=function(x, p) stats::dexp(x, p[1], log=TRUE)),
  dgamma=list(p=c(shape=2.4, rate=1.6), lo=0, hi=Inf,
              ref=function(x, p) stats::dgamma(x, p[1], rate=p[2], log=TRUE)),
  invGamma=list(p=c(alpha=3.2, beta=1.5), lo=0, hi=Inf, ref=NULL),
  dweibull=list(p=c(shape=1.8, scale=2.2), lo=0, hi=Inf,
                ref=function(x, p) stats::dweibull(x, p[1], p[2], log=TRUE)),
  frechet=list(p=c(alpha=2.5, sigma=1.4), lo=0, hi=Inf, ref=NULL),
  rayleigh=list(p=c(sigma=1.3), lo=0, hi=Inf, ref=NULL),
  pareto=list(p=c(y_min=1.2, alpha=3.1), lo=1.2, hi=Inf, ref=NULL),
  paretoType2=list(p=c(mu=0.5, lambda=1.4, alpha=3.3), lo=0.5, hi=Inf, ref=NULL),
  dbeta=list(p=c(shape1=2.1, shape2=3.4), lo=0, hi=1,
             ref=function(x, p) stats::dbeta(x, p[1], p[2], log=TRUE)),
  betaProportion=list(p=c(mu=0.4, kappa=7), lo=0, hi=1,
                      ref=function(x, p) stats::dbeta(x, p[1]*p[2], (1 - p[1])*p[2], log=TRUE)),
  dunif=list(p=c(min=-1, max=2), lo=-1, hi=2,
             ref=function(x, p) stats::dunif(x, p[1], p[2], log=TRUE)))

## The templates are rxode2 expressions; evaluate them in an environment that
## gives the rxode2 spellings their R meaning.  llik*() is rxode2ll's, reached
## through the same arguments in the same order.
.edPeerEnv <- function() {
  list2env(list(
    lgammafn=base::lgamma,
    llikNorm=function(x, m, s) stats::dnorm(x, m, s, log=TRUE),
    llikT=function(x, df, m, s) stats::dt((x - m)/s, df, log=TRUE) - log(s),
    llikCauchy=function(x, l, s) stats::dcauchy(x, l, s, log=TRUE),
    llikChisq=function(x, df) stats::dchisq(x, df, log=TRUE),
    llikExp=function(x, r) stats::dexp(x, r, log=TRUE),
    llikGamma=function(x, sh, rt) stats::dgamma(x, sh, rate=rt, log=TRUE),
    llikWeibull=function(x, sh, sc) stats::dweibull(x, sh, sc, log=TRUE),
    llikBeta=function(x, a, b) stats::dbeta(x, a, b, log=TRUE),
    llikUnif=function(x, a, b) stats::dunif(x, a, b, log=TRUE)),
    parent=baseenv())
}

test_that("every peer log density is the right density, normalized", {
  .defs <- strsplit(nlmixr2est:::.etaDistPeerDefs, "|", fixed=TRUE)
  .env <- .edPeerEnv()
  for (.d in .defs) {
    .nm <- .d[[1]]
    .sp <- .edPeerSpec[[.nm]]
    expect_false(is.null(.sp), label=paste0("test spec for '", .nm, "'"))
    if (is.null(.sp)) next
    .tm <- .d[[2]]
    for (.i in seq_along(.sp$p)) {
      .tm <- gsub(paste0("{", names(.sp$p)[.i], "}"), paste0("(", .sp$p[.i], ")"),
                  .tm, fixed=TRUE)
    }
    .tm <- gsub("{x}", "x", .tm, fixed=TRUE)
    expect_false(grepl("{", .tm, fixed=TRUE),
                 label=paste0("'", .nm, "' template fully substituted"))
    .f <- eval(parse(text=paste0("function(x) { ", .tm, " }")), envir=.env)
    if (!is.null(.sp$ref)) {
      .xs <- if (is.finite(.sp$hi)) {
        seq(.sp$lo + 1e-3, .sp$hi - 1e-3, length.out=7L)
      } else if (.sp$lo == -Inf) {
        seq(-2, 2, length.out=7L)
      } else {
        .sp$lo + c(0.05, 0.3, 0.8, 1.5, 2.5, 4, 6)
      }
      expect_equal(.f(.xs), .sp$ref(.xs, .sp$p), tolerance=1e-10,
                   label=paste0("'", .nm, "' vs R"))
    }
    .i1 <- stats::integrate(function(x) exp(.f(x)), .sp$lo, .sp$hi,
                            rel.tol=1e-8)$value
    expect_equal(.i1, 1, tolerance=1e-6,
                 label=paste0("'", .nm, "' integrates to 1"))
  }
})

test_that("the peer table declines a family it has no density for", {
  skip_if_not(requireNamespace("lotri", quietly=TRUE) &&
              !is.null(utils::getFromNamespace("lotriEtaDists", "lotri")))
  .tab <- nlmixr2est:::.etaDistPeerTable()
  ## stdNormal has no parameters, so it has no M-step and no density here
  expect_true(is.na(.tab$logDensity[.tab$name == "stdNormal"]))
  ## everything else lotri can declare has one
  expect_equal(sum(is.na(.tab$logDensity)), 1L)
})

## The peer's emitted derivatives, checked against central differences on the
## expression the assembler actually produced for Bauer's two-gamma model:
##
##   rx_edll_1_ = llikGamma(ETA[1], exp(-THETA[5]), exp(-(THETA[1] + THETA[5])))
##
## Written out rather than re-derived through symengine so the test states what
## the model is expected to contain.  If the assembler's output changes, this
## does not silently follow it -- it keeps checking the chain rule that the
## declaration `dgamma(shape = 1/exp(lclrv), rate = 1/(exp(lclrv)*exp(lclm)))`
## implies, with lclm = THETA_1_ and lclrv = THETA_5_.

test_that("the peer's theta derivatives match central differences", {
  skip_if_not_installed("rxode2ll")
  .ll <- function(x, sh, rt) rxode2ll::llikGamma(x, sh, rt)$fx
  .dS <- function(x, sh, rt) rxode2ll::llikGamma(x, sh, rt)$dShape
  .dR <- function(x, sh, rt) rxode2ll::llikGamma(x, sh, rt)$dRate
  .eta <- 1.7
  .f <- function(th) .ll(.eta, exp(-th[5]), exp(-(th[1] + th[5])))
  .g1 <- function(th) {
    -1*exp(-(th[1] + th[5]))*.dR(.eta, exp(-th[5]), exp(-(th[1] + th[5])))
  }
  .g5 <- function(th) {
    -1*exp(-th[5])*.dS(.eta, exp(-th[5]), exp(-(th[1] + th[5]))) -
      exp(-(th[1] + th[5]))*.dR(.eta, exp(-th[5]), exp(-(th[1] + th[5])))
  }
  .fd <- function(f, th, j, h = 1e-6) {
    .a <- th; .b <- th; .a[j] <- .a[j] + h; .b[j] <- .b[j] - h
    (f(.a) - f(.b))/(2*h)
  }
  ## the starting values of three of the four gamma arms, which span a 22-fold
  ## range of relative variance
  for (.th in list(c(1.5, 1.5, 0.9, 4.2, -3, -3),
                   c(1.9, 1.8, 0.9, 4.2, -0.6, -0.6),
                   c(1.5, 1.5, 0.9, 4.2, 0.7, 0.7))) {
    expect_equal(.g1(.th), .fd(.f, .th, 1L), tolerance = 1e-6)
    expect_equal(.g5(.th), .fd(.f, .th, 5L), tolerance = 1e-6)
  }
})
