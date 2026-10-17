# A saem fit's objective is computed after the fit (saemControl(logLik =
# FALSE)), by Gaussian quadrature on the kernel's event matrix.  For an IOV fit
# the occasion column is a covariate only of the IOV-rewritten model, so that
# matrix has to be the fit's own, not one rebuilt from the data.
.saemIovObjMod <- function() {
  ini({
    tka <- 0.45
    tcl <- 1
    tv <- 3.45
    add.sd <- 0.7
    eta.v ~ 0.1
    iov.cl ~ 0.1 | occ
  })
  model({
    ka <- exp(tka)
    cl <- exp(tcl + iov.cl)
    v <- exp(tv + eta.v)
    linCmt() ~ add(add.sd)
  })
}

.saemNoIovObjMod <- function() {
  ini({
    tka <- 0.45
    tcl <- 1
    tv <- 3.45
    add.sd <- 0.7
    eta.v ~ 0.1
    eta.cl ~ 0.1
  })
  model({
    ka <- exp(tka)
    cl <- exp(tcl + eta.cl)
    v <- exp(tv + eta.v)
    linCmt() ~ add(add.sd)
  })
}

.saemIovObjData <- function() {
  .d <- nlmixr2data::theo_sd
  .d$occ <- ifelse(.d$TIME < 6, 1L, 2L)
  .d
}

# The fit's objective, its type, and whether reading it fell back to FOCEi
.saemIovObj <- function(mod, iovMethod = NULL) {
  .ctl <- saemControl(print = 0, covMethod = "")
  if (!is.null(iovMethod)) {
    .ctl$iovMethod <- iovMethod
  }
  .fit <- suppressWarnings(suppressMessages(nlmixr2(mod, .saemIovObjData(), "saem", control = .ctl)))
  .msg <- capture.output(.objf <- .fit$objf, type = "message")
  list(objf = .objf, ofvType = .fit$ofvType, fellBack = any(grepl("quadrature failed", .msg)))
}

nmTest({
  test_that("a saem IOV fit reports its Gaussian-quadrature objective", {
    for (.m in c("twoLevel", "theta")) {
      .r <- .saemIovObj(.saemIovObjMod, .m)
      expect_identical(.r$ofvType, "gauss3_1.6", label = .m)
      expect_false(.r$fellBack, label = .m)
      expect_equal(.r$objf, 204.6226, tolerance = 1e-6, label = .m)
    }
    # a fit without IOV is unchanged
    .r <- .saemIovObj(.saemNoIovObjMod)
    expect_identical(.r$ofvType, "gauss3_1.6")
    expect_false(.r$fellBack)
    expect_equal(.r$objf, 181.3438, tolerance = 1e-6)
  })
})
