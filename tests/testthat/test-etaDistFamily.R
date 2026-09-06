nmTest({
  # 22 hand-written quantiles and log densities in src/etaDistFam.cpp back the
  # ODE-free distribution M-step.  That is exactly the kind of code that is
  # silently wrong, so every family is pinned against an independent R
  # reference -- the same reason saemFormGTest() exists.
  #
  # The family code IS the row number in lotri::lotriEtaDists(), so this also
  # catches the C++ dispatch drifting from the catalog.
  .ref <- list(
    dnorm             = list(a=c(2,1.5),     q=function(u,a) qnorm(u,a[1],a[2]),               d=function(x,a) dnorm(x,a[1],a[2],TRUE)),
    stdNormal         = list(a=numeric(0),   q=function(u,a) qnorm(u),                          d=function(x,a) dnorm(x,log=TRUE)),
    studentT          = list(a=c(5,1,2),     q=function(u,a) a[2]+a[3]*qt(u,a[1]),              d=function(x,a) dt((x-a[2])/a[3],a[1],log=TRUE)-log(a[3])),
    dcauchy           = list(a=c(1,2),       q=function(u,a) qcauchy(u,a[1],a[2]),              d=function(x,a) dcauchy(x,a[1],a[2],TRUE)),
    doubleExponential = list(a=c(1,2),       q=function(u,a) a[1]-a[2]*sign(u-.5)*log1p(-2*sign(u-.5)*(u-.5)), d=function(x,a) -log(2*a[2])-abs(x-a[1])/a[2]),
    dlogis            = list(a=c(1,2),       q=function(u,a) qlogis(u,a[1],a[2]),               d=function(x,a) dlogis(x,a[1],a[2],TRUE)),
    gumbel            = list(a=c(1,2),       q=function(u,a) a[1]-a[2]*log(-log(u)),            d=function(x,a) {z<-(x-a[1])/a[2]; -log(a[2])-z-exp(-z)}),
    dlnorm            = list(a=c(0.5,0.8),   q=function(u,a) qlnorm(u,a[1],a[2]),               d=function(x,a) dlnorm(x,a[1],a[2],TRUE)),
    dchisq            = list(a=c(4),         q=function(u,a) qchisq(u,a[1]),                    d=function(x,a) dchisq(x,a[1],log=TRUE)),
    invChiSquare      = list(a=c(6),         q=function(u,a) 1/qchisq(1-u,a[1]),                d=function(x,a) dchisq(1/x,a[1],log=TRUE)-2*log(x)),
    scaledInvChiSquare= list(a=c(6,1.3),     q=function(u,a) a[1]*a[2]^2/qchisq(1-u,a[1]),      d=function(x,a) {nu<-a[1];t2<-a[2]^2; (nu/2)*log(nu*t2/2)-lgamma(nu/2)-(1+nu/2)*log(x)-nu*t2/(2*x)}),
    dexp              = list(a=c(0.7),       q=function(u,a) qexp(u,a[1]),                      d=function(x,a) dexp(x,a[1],TRUE)),
    dgamma            = list(a=c(11.6,2.31), q=function(u,a) qgamma(u,a[1],a[2]),               d=function(x,a) dgamma(x,a[1],a[2],log=TRUE)),
    invGamma          = list(a=c(3,2),       q=function(u,a) a[2]/qgamma(1-u,a[1]),             d=function(x,a) a[1]*log(a[2])-lgamma(a[1])-(a[1]+1)*log(x)-a[2]/x),
    dweibull          = list(a=c(2,3),       q=function(u,a) qweibull(u,a[1],a[2]),             d=function(x,a) dweibull(x,a[1],a[2],TRUE)),
    frechet           = list(a=c(3,2),       q=function(u,a) a[2]*(-log(u))^(-1/a[1]),          d=function(x,a) {z<-x/a[2]; log(a[1]/a[2])-(1+a[1])*log(z)-z^(-a[1])}),
    rayleigh          = list(a=c(1.5),       q=function(u,a) a[1]*sqrt(-2*log1p(-u)),           d=function(x,a) log(x)-2*log(a[1])-x^2/(2*a[1]^2)),
    pareto            = list(a=c(1,3),       q=function(u,a) a[1]*(1-u)^(-1/a[2]),              d=function(x,a) log(a[2])+a[2]*log(a[1])-(a[2]+1)*log(x)),
    paretoType2       = list(a=c(0,2,3),     q=function(u,a) a[1]+a[2]*((1-u)^(-1/a[3])-1),     d=function(x,a) {z<-(x-a[1])/a[2]; log(a[3]/a[2])-(a[3]+1)*log1p(z)}),
    dbeta             = list(a=c(2,5),       q=function(u,a) qbeta(u,a[1],a[2]),                d=function(x,a) dbeta(x,a[1],a[2],log=TRUE)),
    betaProportion    = list(a=c(0.3,10),    q=function(u,a) qbeta(u,a[1]*a[2],(1-a[1])*a[2]),  d=function(x,a) dbeta(x,a[1]*a[2],(1-a[1])*a[2],log=TRUE)),
    dunif             = list(a=c(-1,3),      q=function(u,a) qunif(u,a[1],a[2]),                d=function(x,a) dunif(x,a[1],a[2],TRUE)))

  .tab <- lotri::lotriEtaDists()
  .u <- c(0.01, 0.1, 0.25, 0.5, 0.75, 0.9, 0.99)

  test_that("every catalog family is in the C++ dispatch", {
    expect_true(all(.tab$name %in% names(.ref)))
    for (.i in seq_len(nrow(.tab))) {
      expect_false(is.null(rxEtaDistTest_(.i, .u, .ref[[.tab$name[.i]]]$a)),
                   info = .tab$name[.i])
    }
  })

  test_that("the C++ quantile and log density match R for every family", {
    for (.i in seq_len(nrow(.tab))) {
      .nm <- .tab$name[.i]; .r <- .ref[[.nm]]
      .got <- rxEtaDistTest_(.i, .u, .r$a)
      .qr <- .r$q(.u, .r$a)
      expect_equal(.got$q, .qr, tolerance = 1e-8, info = paste(.nm, "quantile"))
      expect_equal(.got$logd, .r$d(.qr, .r$a), tolerance = 1e-8,
                   info = paste(.nm, "log density"))
      expect_equal(.got$narg, length(.r$a), info = paste(.nm, "narg"))
    }
  })

  test_that("the family code is the catalog row, and unknown falls back", {
    for (.i in seq_len(nrow(.tab))) {
      expect_equal(nlmixr2est:::.etaDistFamilyCode(paste0(.tab$name[.i], "(1,2)")),
                   .i, info = .tab$name[.i])
    }
    # anything not in the catalog gets 0, which selects the R fallback
    expect_equal(nlmixr2est:::.etaDistFamilyCode("dnotafamily(1,2)"), 0L)
  })
})
