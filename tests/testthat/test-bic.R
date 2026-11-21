test_that("BIC() links thetas to etas through intermediate assignments", {
  m <- function() {
    ini({
      tka <- 0.45
      tcl <- log(2.7)
      cwt <- 0.75
      tv <- 3.45
      tf <- 0.1
      eta.cl ~ 0.3
      eta.f ~ 0.1
      prop.sd <- 0.1
      eta.r ~ 0.1
      add.sd <- 0.7
    })
    model({
      ka <- exp(tka)
      tvcl <- exp(tcl) * (WT / 70)^cwt
      cl <- tvcl * exp(eta.cl)
      v <- exp(tv)
      # k carries eta.cl through cl but v has no eta: tv stays fixed-class
      k <- cl / v
      ef <- eta.f
      f(depot) <- tf + ef
      rsd <- prop.sd * exp(eta.r)
      d / dt(depot) <- -ka * depot
      d / dt(center) <- ka * depot - k * center
      cp <- center / v
      cp ~ add(add.sd) + prop(rsd)
    })
  }
  ui <- suppressWarnings(rxode2::rxode2(m))
  .l <- .bicThetaLink(ui)
  expect_setequal(.l$structural, c("tcl", "cwt", "tf"))
  expect_setequal(.l$residual, "prop.sd")
  .c <- .bicParClass(ui)
  expect_equal(.c$total, 10L)
  expect_equal(.c$omega, 3L)
  expect_equal(.c$structural, 3L)
  expect_equal(.c$residual, 1L)
})

test_that("BIC() links only the eta's own sub-expression", {
  m <- function() {
    ini({
      tka <- 0.45
      tcl <- log(2.7)
      tv <- 3.45
      eta.cl ~ 0.3
      add.sd <- 0.7
      prop.sd <- 0.1
      eta.r ~ 0.1
    })
    model({
      ka <- exp(tka)
      v <- exp(tv)
      d / dt(depot) <- -ka * depot
      d / dt(center) <- ka * depot - exp(tcl + eta.cl) * center / v
      cp <- center / v
      w <- cp * prop.sd * exp(eta.r)
      cp ~ add(add.sd) + prop(w)
    })
  }
  .l <- .bicThetaLink(suppressWarnings(rxode2::rxode2(m)))
  expect_setequal(.l$structural, "tcl")
  expect_setequal(.l$residual, "prop.sd")
})

test_that("BIC() stops at a linCmt() prediction and follows indirect residual etas", {
  m <- function() {
    ini({
      tka <- 0.45
      tcl <- log(2.7)
      tv <- 3.45
      eta.cl ~ 0.3
      add.sd <- 0.7
      prop.sd <- 0.1
      eta.r ~ 0.1
    })
    model({
      ka <- exp(tka)
      cl <- exp(tcl + eta.cl)
      v <- exp(tv)
      cp <- linCmt()
      rsd0 <- prop.sd * exp(eta.r)
      rsd <- sqrt(add.sd^2 + (rsd0 * cp)^2)
      cp ~ add(rsd)
    })
  }
  .l <- .bicThetaLink(suppressWarnings(rxode2::rxode2(m)))
  expect_setequal(.l$structural, "tcl")
  expect_setequal(.l$residual, "prop.sd")
})

test_that("BIC() treats a ll() endpoint as structural", {
  m <- function() {
    ini({
      t1 <- 0.1
      t2 <- 0.2
      eta1 ~ 0.1
    })
    model({
      lp <- t1 + eta1
      ll(err) ~ DV * lp - log(1 + exp(lp)) + t2 * 0
    })
  }
  .l <- .bicThetaLink(suppressWarnings(rxode2::rxode2(m)))
  expect_setequal(.l$structural, "t1")
  expect_length(.l$residual, 0L)
})

test_that("BIC() follows rescaled etas, sums and algebraic predictions", {
  .link <- function(m) .bicThetaLink(suppressWarnings(rxode2::rxode2(m)))
  # a theta applied after the eta is part of the same parameter
  m1 <- function() {
    ini({
      tcl <- log(2.7)
      cwt <- 0.75
      eta.cl ~ 0.3
      add.sd <- 0.7
    })
    model({
      cl0 <- exp(tcl + eta.cl)
      cl <- cl0 * (WT / 70)^cwt
      cp <- cl * t
      cp ~ add(add.sd)
    })
  }
  .l <- .link(m1)
  expect_setequal(.l$structural, c("tcl", "cwt"))
  expect_length(.l$residual, 0L)
  # an algebraic prediction feeding the error weight stays structural; a
  # random term may use a derived covariate
  m2 <- function() {
    ini({
      temax <- 1
      tec50 <- 0.5
      eta.emax ~ 0.1
      eta.w ~ 0.1
      add.sd <- 0.7
      prop.sd <- 0.1
    })
    model({
      lwt <- log(WT / 70)
      ew <- eta.w * lwt
      emax <- temax + eta.emax
      ec50 <- exp(tec50 + ew)
      pred <- emax * Conc / (ec50 + Conc)
      y <- pred
      rsd <- sqrt(add.sd^2 + (prop.sd * pred)^2)
      y ~ add(rsd)
    })
  }
  .l <- .link(m2)
  expect_setequal(.l$structural, c("temax", "tec50"))
  expect_length(.l$residual, 0L)
  # only the eta-bearing term of a sum is tied to the eta
  m3 <- function() {
    ini({
      tcl <- log(2.7)
      tka <- 0.45
      eta.cl ~ 0.3
      add.sd <- 0.7
      prop.sd <- 0.1
      eta.r ~ 0.1
    })
    model({
      k <- exp(tka) + exp(tcl + eta.cl)
      cp <- k * t
      rsd <- add.sd + prop.sd * exp(eta.r)
      cp ~ prop(rsd)
    })
  }
  .l <- .link(m3)
  expect_setequal(.l$structural, "tcl")
  expect_setequal(.l$residual, "prop.sd")
})

test_that("BIC() keeps separately assigned parameters out of an eta's product", {
  m <- function() {
    ini({
      te0 <- 1
      temax <- 1
      tec50 <- 0.5
      tcl <- log(2.7)
      tv <- 3.45
      eta.e0 ~ 0.1
      eta.cl ~ 0.3
      add.sd <- 0.7
    })
    model({
      e0 <- te0 + eta.e0
      emax <- exp(temax)
      ec50 <- exp(tec50)
      v <- exp(tv)
      k <- exp(tcl + eta.cl) / v / exp(tv)
      effect <- e0 + emax * Conc / (ec50 + Conc) - k * t
      effect ~ add(add.sd)
    })
  }
  .l <- .bicThetaLink(suppressWarnings(rxode2::rxode2(m)))
  expect_setequal(.l$structural, c("te0", "tcl"))
  expect_length(.l$residual, 0L)
})

test_that("BIC() keeps a mu-referenced residual theta residual", {
  m <- function() {
    ini({
      tcl <- log(2.7)
      eta.cl ~ 0.3
      tprop <- log(0.1)
      eta.rsd ~ 0.1
    })
    model({
      cl <- exp(tcl + eta.cl)
      cp <- cl * t
      rsd <- exp(tprop + eta.rsd)
      cp ~ prop(rsd)
    })
  }
  .c <- .bicParClass(suppressWarnings(rxode2::rxode2(m)))
  expect_equal(.c$structural, 1L)
  expect_equal(.c$residual, 1L)
})

test_that("BIC() does not count same() copy rows", {
  m <- function() {
    ini({
      tcl <- log(2.7)
      eta.cl ~ 0.3
      eta.cl2 ~ 0.3
      add.sd <- 0.7
    })
    model({
      cp <- exp(tcl + eta.cl) * exp(eta.cl2) * t
      cp ~ add(add.sd)
    })
  }
  ui <- suppressWarnings(rxode2::rxode2(m))
  .idf <- ui$iniDf
  .idf$condition[.idf$name == "eta.cl2"] <- "id:same:eta.cl"
  .u <- list(
    iniDf = .idf,
    lstExpr = ui$lstExpr,
    state = ui$state,
    muRefDataFrame = ui$muRefDataFrame
  )
  .c <- .bicParClass(.u)
  expect_equal(.c$total, 3L)
  expect_equal(.c$omega, 1L)
})

test_that("BIC() does not count fixed parameters", {
  m <- function() {
    ini({
      tka <- 0.45
      tcl <- fix(1)
      tv <- 3.45
      eta.ka ~ fix(0.6)
      eta.cl ~ 0.3
      eta.v ~ 0.1
      add.sd <- 0.7
    })
    model({
      ka <- exp(tka + eta.ka)
      cl <- exp(tcl + eta.cl)
      v <- exp(tv + eta.v)
      linCmt() ~ add(add.sd)
    })
  }
  .c <- .bicParClass(suppressWarnings(rxode2::rxode2(m)))
  expect_equal(.c$total, 5L)
  expect_equal(.c$omega, 2L)
  expect_equal(.c$structural, 2L)
  expect_equal(.c$residual, 0L)
})

nmTest({
  test_that("BIC() variants match hand-computed values", {
    one.cmt <- function() {
      ini({
        tka <- 0.45
        tcl <- log(2.7)
        tv <- 3.45
        eta.ka ~ 0.6
        eta.cl ~ 0.3
        add.sd <- 0.7
      })
      model({
        ka <- exp(tka + eta.ka)
        cl <- exp(tcl + eta.cl)
        v <- exp(tv)
        linCmt() ~ add(add.sd)
      })
    }
    fit <- .nlmixr(one.cmt, theo_sd, est = "focei", control = foceiControlFast)
    .m2ll <- -2 * as.numeric(logLik(fit))
    .lo <- log(132)
    .ls <- log(12)
    expect_equal(BIC(fit), fit$BIC)
    expect_equal(BIC(fit, type = "fixed"), fit$BIC)
    expect_equal(BIC(fit, type = "fixed"), .m2ll + 6 * .lo)
    expect_equal(BIC(fit, type = "iiv"), .m2ll + 4 * .lo + 2 * .ls)
    expect_equal(BIC(fit, type = "mixed"), .m2ll + 2 * .lo + 4 * .ls)
    expect_equal(BIC(fit, type = "random"), .m2ll + 2 * .lo + 4 * .ls)
    expect_equal(attr(.bicValue(fit, "mixed"), "dfSub"), 4)
    expect_error(BIC(fit, type = "bad"))
    expect_equal(BIC(fit, "mixed"), BIC(fit, type = "mixed"))
    fit2 <- fit
    .two <- BIC(fit, fit2, type = "mixed")
    expect_equal(names(.two), c("df", "BIC"))
    expect_equal(row.names(.two), c("fit", "fit2"))
    expect_equal(BIC(fit, fit2, "mixed"), .two)
    expect_equal(.two$BIC, rep(.m2ll + 2 * .lo + 4 * .ls, 2))
    expect_equal(.two$df, c(6, 6))
  })

  test_that("BIC() counts inter-occasion kappas as random-effect parameters", {
    one.cmt <- function() {
      ini({
        tka <- 0.45
        tcl <- log(2.7)
        tv <- 3.45
        eta.cl ~ 0.3
        iov.cl ~ 0.1 | occ
        add.sd <- 0.7
      })
      model({
        ka <- exp(tka)
        cl <- exp(tcl + eta.cl + iov.cl)
        v <- exp(tv)
        linCmt() ~ add(add.sd)
      })
    }
    d <- nlmixr2data::theo_md
    d$occ <- 1
    d$occ[d$TIME >= 144] <- 2
    fit <- .nlmixr(one.cmt, d, est = "focei", control = foceiControlFast)
    .m2ll <- -2 * as.numeric(logLik(fit))
    .lo <- log(fit$nobs)
    .ls <- log(12)
    expect_equal(attr(logLik(fit), "df"), 6)
    expect_equal(BIC(fit, type = "fixed"), BIC(fit))
    expect_equal(BIC(fit, type = "iiv"), .m2ll + 4 * .lo + 2 * .ls)
    expect_equal(BIC(fit, type = "mixed"), .m2ll + 3 * .lo + 3 * .ls)
  })

  test_that("BIC() charges a residual-error theta log(nsub) only for type='random'", {
    one.cmt <- function() {
      ini({
        tka <- 0.45
        tcl <- log(2.7)
        tv <- 3.45
        eta.cl ~ 0.3
        add.sd <- 0.7
        eta.r ~ 0.1
      })
      model({
        ka <- exp(tka)
        cl <- exp(tcl + eta.cl)
        v <- exp(tv)
        rsd <- add.sd * exp(eta.r)
        linCmt() ~ add(rsd)
      })
    }
    fit <- .nlmixr(one.cmt, theo_sd, est = "focei", control = foceiControlFast)
    .m2ll <- -2 * as.numeric(logLik(fit))
    .lo <- log(132)
    .ls <- log(12)
    expect_equal(BIC(fit, type = "mixed"), .m2ll + 3 * .lo + 3 * .ls)
    expect_equal(BIC(fit, type = "random"), .m2ll + 2 * .lo + 4 * .ls)
  })
})
