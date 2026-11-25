# nlmixr2ObjGrad(): the fit's objective and analytic outer gradient at
# user-supplied parameters, checked against central differences of that objective.

nmTest({
  .ogModel <- function() {
    ini({
      tka <- 0.45
      tcl <- 1
      tv <- 3.45
      eta.ka ~ 0.6
      eta.cl + eta.v ~ c(0.3, 0.05, 0.1)
      add.sd <- 0.7
    })
    model({
      ka <- exp(tka + eta.ka)
      cl <- exp(tcl + eta.cl)
      v <- exp(tv + eta.v)
      d / dt(depot) <- -ka * depot
      d / dt(center) <- ka * depot - cl / v * center
      cp <- center / v
      cp ~ add(add.sd)
    })
  }

  .ogCtl <- function(fast = TRUE, ...) {
    foceiControl(print = 0L, covMethod = "", calcTables = FALSE, fast = fast, ...)
  }

  # central difference of the objective in theta and in the omega elements
  .ogFd <- function(fit, theta, omega, h = 1e-4) {
    .obj <- function(th, om) nlmixr2ObjGrad(fit, theta = th, omega = om)$objective
    .gt <- vapply(
      names(theta),
      function(n) {
        tp <- theta
        tm <- theta
        tp[n] <- tp[n] + h
        tm[n] <- tm[n] - h
        (.obj(tp, omega) - .obj(tm, omega)) / (2 * h)
      },
      numeric(1)
    )
    .ini <- fit$ui$iniDf
    .ini <- .ini[!is.na(.ini$neta1) & !.ini$fix, , drop = FALSE]
    .go <- vapply(
      seq_len(nrow(.ini)),
      function(k) {
        i <- .ini$neta1[k]
        j <- .ini$neta2[k]
        op <- omega
        om <- omega
        op[i, j] <- op[j, i] <- omega[i, j] + h
        om[i, j] <- om[j, i] <- omega[i, j] - h
        (.obj(theta, op) - .obj(theta, om)) / (2 * h)
      },
      numeric(1)
    )
    c(.gt, stats::setNames(.go, .ini$name))
  }

  test_that("nlmixr2ObjGrad() validates its arguments", {
    expect_error(nlmixr2ObjGrad(list()), "nlmixr2 fit")
  })

  test_that("nlmixr2ObjGrad() reproduces the fit and matches central differences", {
    skip_on_cran()
    skip_if_not_installed("nlmixr2data")
    .th <- rxode2::rxCores()
    on.exit(rxode2::setRxThreads(.th), add = TRUE)
    rxode2::setRxThreads(1L)
    fit <- suppressMessages(nlmixr2(.ogModel, nlmixr2data::theo_sd, "focei", .ogCtl()))

    # at the fit: the fit's own objective, and a ~zero gradient
    og <- nlmixr2ObjGrad(fit, scale = "transformed")
    expect_equal(og$objective, fit$objective, tolerance = 1e-8)
    expect_gt(attr(og, "nAnalyticGradDirect"), 0L)
    expect_true(all(abs(og$gradient) < 0.5))
    expect_equal(og$par[names(fit$theta)], fit$theta)

    # off the fit, untransformed scale: thetas and omega elements against central
    # differences
    th <- fit$theta
    th["tka"] <- th["tka"] + 0.1
    th["add.sd"] <- th["add.sd"] * 1.1
    om <- fit$omega * 1.2
    og <- nlmixr2ObjGrad(fit, theta = th, omega = om)
    expect_gt(attr(og, "nAnalyticGradDirect"), 0L)
    expect_equal(og$scale, "untransformed")
    expect_setequal(names(og$gradient), c(names(th), "eta.ka", "eta.cl", "(eta.cl,eta.v)", "eta.v"))
    fd <- .ogFd(fit, th, om)
    expect_equal(og$gradient[names(fd)], fd, tolerance = 2e-3)
    expect_equal(og$par[["(eta.cl,eta.v)"]], om["eta.cl", "eta.v"])

    # a named subset of theta leaves the rest at the fit's values
    og1 <- nlmixr2ObjGrad(fit, theta = th["tka"])
    th1 <- fit$theta
    th1["tka"] <- th["tka"]
    expect_equal(og1$objective, nlmixr2ObjGrad(fit, theta = unname(th1))$objective)

    # transformed scale: the omega gradient is with respect to the estimation-scale
    # (Cholesky) parameters, checked by differencing through them
    ogt <- nlmixr2ObjGrad(fit, theta = th, omega = om, scale = "transformed")
    expect_equal(ogt$objective, og$objective)
    expect_equal(ogt$gradient[names(th)], og$gradient[names(th)])
    .inv <- rxode2::rxSymInvCholCreate(mat = om, diag.xform = "sqrt")
    .ch0 <- .inv$theta
    .omAt <- function(ch) {
      .inv$theta <- ch
      .m <- .inv$omega
      dimnames(.m) <- dimnames(om)
      .m
    }
    .k <- grep("^om[.]chol[.]", names(ogt$gradient))
    expect_length(.k, 4L)
    h <- 1e-4
    fdc <- vapply(
      seq_along(.ch0),
      function(k) {
        cp <- .ch0
        cm <- .ch0
        cp[k] <- cp[k] + h
        cm[k] <- cm[k] - h
        (nlmixr2ObjGrad(fit, theta = th, omega = .omAt(cp))$objective -
          nlmixr2ObjGrad(fit, theta = th, omega = .omAt(cm))$objective) /
          (2 * h)
      },
      numeric(1)
    )
    expect_equal(unname(ogt$gradient[.k]), fdc, tolerance = 2e-3)
    expect_equal(unname(ogt$par[.k]), .ch0)

    expect_error(nlmixr2ObjGrad(fit, theta = c(nope = 1)), "unknown")
    expect_error(nlmixr2ObjGrad(fit, theta = 1:2), "length")
    expect_error(nlmixr2ObjGrad(fit, omega = diag(2)), "dimensions")
  })

  test_that("nlmixr2ObjGrad() works for a fast=FALSE fit and under est='foce'", {
    skip_on_cran()
    skip_if_not_installed("nlmixr2data")
    # away from the optimum the reference needs a tight inner problem: every
    # differenced objective re-optimizes the etas
    fit <- suppressMessages(nlmixr2(
      .ogModel,
      nlmixr2data::theo_sd,
      "foce",
      .ogCtl(
        fast = FALSE,
        maxOuterIterations = 0L,
        sigdig = 4,
        maxInnerIterations = 300L,
        epsilon = 1e-10,
        trustFterm = 1e-12,
        trustMterm = 1e-12
      )
    ))
    th <- fit$theta
    og <- nlmixr2ObjGrad(fit)
    expect_gt(attr(og, "nAnalyticGradDirect"), 0L)
    expect_equal(og$objective, fit$objective, tolerance = 1e-8)
    fd <- .ogFd(fit, th, fit$omega, h = 1e-3)
    expect_equal(og$gradient[names(fd)], fd, tolerance = 1e-2)
  })

  test_that("nlmixr2ObjGrad() is right on a multi-endpoint model", {
    skip_on_cran()
    me <- function() {
      ini({
        tka <- 0.5
        tcl <- -2
        tv <- 2
        tbase <- 1
        tkout <- -1
        add.pk <- 1
        add.pd <- 0.5
        eta.cl ~ 0.1
      })
      model({
        ka <- exp(tka)
        cl <- exp(tcl + eta.cl)
        v <- exp(tv)
        d / dt(depot) <- -ka * depot
        d / dt(center) <- ka * depot - cl / v * center
        cp <- center / v
        pd <- exp(tbase) * cp / (exp(tkout) + cp)
        cp ~ add(add.pk) | cp
        pd ~ add(add.pd) | pd
      })
    }
    withr::local_seed(7)
    mk <- function(i, mu) {
      o <- rbind(
        data.frame(ID = i, TIME = c(.5, 1, 2, 4, 8), EVID = 0, AMT = 0, DV = abs(stats::rnorm(5, mu, 1)), CMT = "cp"),
        data.frame(
          ID = i,
          TIME = c(.75, 1.5, 3, 6, 10),
          EVID = 0,
          AMT = 0,
          DV = abs(stats::rnorm(5, 2, .3)),
          CMT = "pd"
        )
      )
      rbind(data.frame(ID = i, TIME = 0, EVID = 101, AMT = 100, DV = NA, CMT = "cp"), o[order(o$TIME), ])
    }
    d <- do.call(rbind, lapply(1:4, function(i) mk(i, 3 + i)))
    fit <- suppressWarnings(suppressMessages(
      nlmixr2(me, d, "focei", .ogCtl(maxOuterIterations = 0L))
    ))
    og <- nlmixr2ObjGrad(fit)
    expect_gt(attr(og, "nAnalyticGradDirect"), 0L)
    expect_equal(og$objective, fit$objective, tolerance = 1e-8)
    fd <- .ogFd(fit, fit$theta, fit$omega)
    expect_equal(og$gradient[names(fd)], fd, tolerance = 2e-3)
  })
})
