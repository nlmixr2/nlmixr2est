# Analytic outer gradient (foceiControl(fast=TRUE)) for inter-occasion variability
# models (#1377): every component -- thetas, sigma, between-subject omega and the
# occasion-level variability -- agrees with central differences of the objective,
# under both IOV expansions, and a fast fit lands on the finite-difference fit.

nmTest({
  .iovGradMod <- function() {
    ini({
      tka <- 0.3
      tcl <- 1.2
      tv <- 3.2
      add.sd <- 0.9
      eta.cl ~ 0.2
      eta.v ~ 0.1
      iov.cl ~ 0.05 | occ
    })
    model({
      ka <- exp(tka)
      cl <- exp(tcl + eta.cl + iov.cl)
      v <- exp(tv + eta.v)
      d / dt(depot) <- -ka * depot
      d / dt(center) <- ka * depot - cl / v * center
      cp <- center / v
      cp ~ add(add.sd)
    })
  }

  .iovGradCtl <- function(fast, iovMethod = "auto") {
    foceiControl(
      print = 0L,
      covMethod = "",
      fast = fast,
      sigdig = 8,
      iovMethod = iovMethod,
      rxControl = rxode2::rxControl(atol = 1e-12, rtol = 1e-12),
      maxOuterIterations = 0L,
      maxInnerIterations = 300L,
      calcTables = FALSE
    )
  }

  # central differences of the objective in each ini() value of `ui`; `omega`
  # names are variances (or a "(a,b)" covariance) of an omega block
  .iovCentral <- function(ui, d, ctl, thetas, omega = character(0), h = 1e-4) {
    .obj <- function(u) suppressMessages(suppressWarnings(nlmixr2(u, d, "focei", ctl)))$objf
    .ini <- ui$iniDf
    .at <- function(nm, v) {
      if (nm %in% thetas) {
        return(do.call(rxode2::ini, c(list(ui), stats::setNames(list(v), nm))))
      }
      # rebuild the element's whole block from iniDf, then restate it
      .w <- which(.ini$name == nm)
      .lvl <- .ini$condition[.w]
      .blk <- .ini[is.na(.ini$ntheta) & .ini$condition == .lvl, , drop = FALSE]
      # the block holds the element's etas and every eta correlated with them
      .own <- c(.ini$neta1[.w], .ini$neta2[.w])
      .off <- .blk[.blk$neta1 != .blk$neta2 & (.blk$neta1 %in% .own | .blk$neta2 %in% .own), , drop = FALSE]
      .etas <- unique(c(.own, .off$neta1, .off$neta2))
      .dg <- .blk[.blk$neta1 == .blk$neta2 & .blk$neta1 %in% .etas, , drop = FALSE]
      .dg <- .dg[order(.dg$neta1), , drop = FALSE]
      .m <- diag(.dg$est, nrow(.dg))
      for (.r in which(.blk$neta1 != .blk$neta2)) {
        .a <- match(.blk$neta1[.r], .dg$neta1)
        .b <- match(.blk$neta2[.r], .dg$neta1)
        if (!anyNA(c(.a, .b))) .m[.a, .b] <- .m[.b, .a] <- .blk$est[.r]
      }
      .a <- match(.ini$neta1[.w], .dg$neta1)
      .b <- match(.ini$neta2[.w], .dg$neta1)
      .m[.a, .b] <- .m[.b, .a] <- v
      .lhs <- str2lang(paste(.dg$name, collapse = "+"))
      # column-major upper triangle is lotri's row-major lower triangle
      .rhs <- as.call(c(list(as.name("c")), as.list(.m[upper.tri(.m, diag = TRUE)])))
      if (.lvl != "id") {
        .rhs <- bquote(.(.rhs) | .(as.name(.lvl)))
      }
      eval(bquote(rxode2::ini(ui, .(.lhs) ~ .(.rhs))))
    }
    vapply(
      c(thetas, omega),
      function(nm) {
        .v <- .ini$est[.ini$name == nm]
        (.obj(.at(nm, .v + h)) - .obj(.at(nm, .v - h))) / (2 * h)
      },
      numeric(1)
    )
  }

  # d(om.chol)/d(omega element) for one block, numerically through the same
  # symbolic inverse-Cholesky the fit uses (diagXform="sqrt", the default)
  .iovCholJac <- function(m, h = 1e-7) {
    .th <- function(x) rxode2::rxSymInvCholCreate(mat = x, diag.xform = "sqrt")$theta
    .idx <- which(lower.tri(m, diag = TRUE), arr.ind = TRUE)
    .idx <- .idx[order(.idx[, 1], .idx[, 2]), , drop = FALSE]
    vapply(
      seq_len(nrow(.idx)),
      function(k) {
        .up <- .dn <- m
        .i <- .idx[k, 1]
        .j <- .idx[k, 2]
        .up[.i, .j] <- .up[.j, .i] <- m[.i, .j] + h
        .dn[.i, .j] <- .dn[.j, .i] <- m[.i, .j] - h
        (.th(.up) - .th(.dn)) / (2 * h)
      },
      numeric(length(.th(m)))
    )
  }

  test_that("IOV analytic gradient matches central differences (iovMethod theta)", {
    skip_on_cran()
    skip_if_not_installed("nlmixr2data")
    d <- .iovFitData()
    ph <- suppressMessages(nlmixr2(.iovGradMod, d, "focei", .iovGradCtl(TRUE, "theta")))
    expect_gt(ph$env$nAnalyticGradDirect, 0L)
    g <- .foceiGradDirect(ph)
    expect_identical(names(g), c("tka", "tcl", "tv", "add.sd", "iov.cl", "om.chol.1", "om.chol.2"))
    ui <- ph$finalUi
    fd <- .gradRef("focei-iov-theta", function() {
      .iovCentral(ui, d, .iovGradCtl(FALSE, "theta"), c("tka", "tcl", "tv", "add.sd"), c("eta.cl", "eta.v", "iov.cl"))
    })
    expect_equal(unname(g[c("tka", "tcl", "tv", "add.sd")]), unname(fd[1:4]), tolerance = 1e-3)
    # a diagonal om.chol is var^(-1/4), so d var/d om.chol = -4 var^(5/4)
    .var <- c(eta.cl = 0.2, eta.v = 0.1)
    expect_equal(unname(g[c("om.chol.1", "om.chol.2")]), unname(fd[names(.var)] * -4 * .var^1.25), tolerance = 1e-3)
    # the occasion variance rides on its magnitude theta, an sd
    expect_equal(unname(g["iov.cl"]), unname(fd["iov.cl"] * 2 * sqrt(0.05)), tolerance = 1e-3)
  })

  test_that("correlated IOV analytic gradient matches central differences (iovMethod omega)", {
    skip_on_cran()
    skip_if_not_installed("nlmixr2data")
    mod <- function() {
      ini({
        tka <- 0.3
        tcl <- 1.2
        tv <- 3.2
        add.sd <- 0.9
        eta.cl ~ 0.2
        eta.v ~ 0.1
        iov.cl + iov.v ~ c(0.05, 0.01, 0.04) | occ
      })
      model({
        ka <- exp(tka)
        cl <- exp(tcl + eta.cl + iov.cl)
        v <- exp(tv + eta.v + iov.v)
        d / dt(depot) <- -ka * depot
        d / dt(center) <- ka * depot - cl / v * center
        cp <- center / v
        cp ~ add(add.sd)
      })
    }
    d <- .iovFitData()
    ph <- suppressMessages(nlmixr2(mod, d, "focei", .iovGradCtl(TRUE)))
    expect_gt(ph$env$nAnalyticGradDirect, 0L)
    g <- .foceiGradDirect(ph)
    # two subject-level variances, then ONE occasion block shared by every occasion
    expect_identical(names(g), c("tka", "tcl", "tv", "add.sd", paste0("om.chol.", 1:5)))
    ui <- ph$finalUi
    .om <- c("eta.cl", "eta.v", "iov.cl", "(iov.cl,iov.v)", "iov.v")
    fd <- .gradRef("focei-iov-omega-cor", function() {
      .iovCentral(ui, d, .iovGradCtl(FALSE), c("tka", "tcl", "tv", "add.sd"), .om)
    })
    expect_equal(unname(g[c("tka", "tcl", "tv", "add.sd")]), unname(fd[1:4]), tolerance = 1e-3)
    .var <- c(eta.cl = 0.2, eta.v = 0.1)
    expect_equal(unname(g[c("om.chol.1", "om.chol.2")]), unname(fd[names(.var)] * -4 * .var^1.25), tolerance = 1e-3)
    # occasion block: dObj/d(element) = sum_k dObj/d(om.chol_k) d(om.chol_k)/d(element)
    .jac <- .iovCholJac(matrix(c(0.05, 0.01, 0.01, 0.04), 2L))
    expect_equal(
      as.numeric(crossprod(.jac, g[paste0("om.chol.", 3:5)])),
      unname(fd[.om[3:5]]),
      tolerance = 1e-3
    )
  })

  test_that("IOV with an ll() endpoint takes the analytic gradient", {
    skip_on_cran()
    skip_if_not_installed("nlmixr2data")
    mod <- function() {
      ini({
        tka <- 0.3
        tcl <- 1.2
        tv <- 3.2
        add.sd <- 0.9
        eta.cl ~ 0.2
        eta.v ~ 0.1
        iov.cl ~ 0.05 | occ
      })
      model({
        ka <- exp(tka)
        cl <- exp(tcl + eta.cl + iov.cl)
        v <- exp(tv + eta.v)
        d / dt(depot) <- -ka * depot
        d / dt(center) <- ka * depot - cl / v * center
        cp <- center / v
        ll(err) ~ -0.5 * log(2 * pi) - log(add.sd) - 0.5 * ((DV - cp) / add.sd)^2
      })
    }
    d <- .iovFitData()
    ph <- suppressMessages(nlmixr2(mod, d, "focei", .iovGradCtl(TRUE, "theta")))
    expect_gt(ph$env$nAnalyticGradDirect, 0L)
    g <- .foceiGradDirect(ph)
    ui <- ph$finalUi
    fd <- .gradRef("focei-iov-ll", function() {
      # h = 1e-4 sits in the inner solve's noise for this log-density objective
      .iovCentral(ui, d, .iovGradCtl(FALSE, "theta"), c("tka", "tcl", "tv", "add.sd"), "iov.cl", h = 3e-4)
    })
    expect_equal(unname(g[c("tka", "tcl", "tv", "add.sd")]), unname(fd[1:4]), tolerance = 1e-3)
    expect_equal(unname(g["iov.cl"]), unname(fd["iov.cl"] * 2 * sqrt(0.05)), tolerance = 1e-3)
  })

  test_that("a fast IOV fit lands on the finite-difference fit", {
    skip_on_cran()
    skip_if_not_installed("nlmixr2data")
    d <- .iovFitData()
    .ctl <- function(fast) foceiControl(print = 0L, covMethod = "", fast = fast, calcTables = FALSE)
    fF <- .numRef("fit-fd-iov-ode", function() {
      .f <- suppressMessages(nlmixr2(.iovGradMod, d, "focei", .ctl(FALSE)))
      list(objf = .f$objf, theta = unname(.f$theta))
    })
    fT <- suppressMessages(nlmixr2(.iovGradMod, d, "focei", .ctl(TRUE)))
    expect_gt(fT$env$nAnalyticGradDirect, 0L)
    expect_equal(as.numeric(fT$objf), as.numeric(fF$objf), tolerance = 1e-3)
    expect_equal(unname(fT$theta), fF$theta, tolerance = 1e-2)
  })
})
