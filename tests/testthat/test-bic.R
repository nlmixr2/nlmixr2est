test_that("bic() links thetas to etas through intermediate assignments", {
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

test_that("bic() does not count fixed parameters", {
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
  test_that("bic() variants match hand-computed values", {
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
    expect_equal(as.numeric(bic(fit)), BIC(fit))
    expect_equal(as.numeric(bic(fit, "fixed")), fit$BIC)
    expect_equal(as.numeric(bic(fit, "fixed")), .m2ll + 6 * .lo)
    expect_equal(as.numeric(bic(fit, "iiv")), .m2ll + 4 * .lo + 2 * .ls)
    expect_equal(as.numeric(bic(fit, "mixed")), .m2ll + 2 * .lo + 4 * .ls)
    expect_equal(as.numeric(bic(fit, "random")), .m2ll + 2 * .lo + 4 * .ls)
    expect_equal(attr(bic(fit, "mixed"), "dfSub"), 4)
    expect_error(bic(fit, "bad"))
  })

  test_that("bic() counts inter-occasion kappas as random-effect parameters", {
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
    expect_equal(as.numeric(bic(fit, "fixed")), BIC(fit))
    expect_equal(as.numeric(bic(fit, "iiv")), .m2ll + 4 * .lo + 2 * .ls)
    expect_equal(as.numeric(bic(fit, "mixed")), .m2ll + 3 * .lo + 3 * .ls)
  })

  test_that("bic() charges a residual-error theta log(nsub) only for type='random'", {
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
    expect_equal(as.numeric(bic(fit, "mixed")), .m2ll + 3 * .lo + 3 * .ls)
    expect_equal(as.numeric(bic(fit, "random")), .m2ll + 2 * .lo + 4 * .ls)
  })
})
