.saemEtaMatTestCfg <- function() {
  .n <- 3L
  .nmc <- 2L
  # phi columns: 1 and 3 carry etas (phi1), 2 is population-only (phi0)
  .phiM <- matrix(seq_len(.n * .nmc * 3) / 10, .n * .nmc, 3)
  list(
    N = .n,
    nmc = .nmc,
    i1 = c(0L, 2L),
    i0 = 1L,
    mprior_phi1 = cbind(rep(1, .n), rep(2, .n)),
    phiM = .phiM,
    ue = matrix(1, .n * .nmc, 3),
    statphi11 = NULL,
    statphi01 = NULL,
    statphi12 = NULL,
    statphi02 = NULL
  )
}

.saemEtaMatTestUi <- function(mixProbs = NULL) {
  list(
    iniDf = data.frame(
      name = c("tka", "tv", "eta.ka", "eta.v"),
      neta1 = c(NA, NA, 1, 2),
      neta2 = c(NA, NA, 1, 2)
    ),
    # eta.ka -> phi1 column 2, eta.v -> phi1 column 1
    saemOmegaTrans = c(2L, 1L),
    mixProbs = mixProbs
  )
}

test_that("saemControl(etaMat=) normalises its argument", {
  expect_identical(saemControl(), saemControl(etaMat = NULL))
  expect_null(saemControl()$etaMat)
  .m <- matrix(1:6, 3, 2, dimnames = list(NULL, c("a", "b")))
  expect_identical(saemControl(etaMat = .m)$etaMat, matrix(as.double(1:6), 3, 2, dimnames = list(NULL, c("a", "b"))))
  .df <- data.frame(ID = c(10, 20), a = c(0.1, 0.2))
  expect_identical(rownames(saemControl(etaMat = .df)$etaMat), c("10", "20"))
  expect_identical(colnames(saemControl(etaMat = .df)$etaMat), "a")
  # a stored control re-run through saemControl() keeps its IDs
  .ctl <- saemControl(etaMat = .df)
  expect_identical(do.call(saemControl, list(etaMat = .ctl$etaMat))$etaMat, .ctl$etaMat)
  expect_error(saemControl(etaMat = matrix(NA_real_, 2, 2)), "etaMat")
  expect_error(saemControl(etaMat = matrix(Inf, 2, 2)), "finite")
  expect_error(saemControl(etaMat = data.frame(ID = c(1, 1), a = 1:2)), "duplicated ID")
  expect_error(saemControl(etaMat = "a"), "etaMat")
})

test_that(".saemEtaMatCfg() starts every chain at the seed", {
  .cfg <- .saemEtaMatTestCfg()
  .ui <- .saemEtaMatTestUi()
  expect_identical(.saemEtaMatCfg(.cfg, NULL, .ui, 1:3), .cfg)
  .eta <- cbind(eta.ka = c(0.1, 0.2, 0.3), eta.v = c(-0.1, -0.2, -0.3))
  .ret <- .saemEtaMatCfg(.cfg, .eta, .ui, 1:3)
  .rows <- rep(1:3, 2)
  expect_equal(.ret$phiM[, 1], 1 + .eta[.rows, "eta.v"])
  expect_equal(.ret$phiM[, 3], 2 + .eta[.rows, "eta.ka"])
  # the population-only column is untouched
  expect_equal(.ret$phiM[, 2], .cfg$phiM[, 2])
  expect_equal(.ret$statphi11, cbind(1 + .eta[, "eta.v"], 2 + .eta[, "eta.ka"]))
  expect_equal(.ret$statphi12, crossprod(.ret$phiM[, c(1, 3)]))
  # unnamed columns follow the model's eta order; IDs in row names are matched
  .ret2 <- .saemEtaMatCfg(.cfg, unname(.eta), .ui, 1:3)
  expect_equal(.ret2$phiM, .ret$phiM)
  .idEta <- .eta[3:1, ]
  rownames(.idEta) <- c("c", "b", "a")
  expect_equal(.saemEtaMatCfg(.cfg, .idEta, .ui, c("a", "b", "c"))$phiM, .ret$phiM)
  .ret3 <- .saemEtaMatCfg(.cfg, .eta[3:1, ], .ui, 1:3)
  expect_equal(.ret3$phiM[1, 1], 1 + .eta[[3, "eta.v"]])
  # columns in another order are matched by name
  expect_equal(.saemEtaMatCfg(.cfg, .eta[, 2:1], .ui, 1:3)$phiM, .ret$phiM)
})

test_that(".saemEtaMatCfg() keeps uninformative and temporary etas at the prior", {
  .cfg <- .saemEtaMatTestCfg()
  .cfg$ue[, 3] <- 0
  .ui <- .saemEtaMatTestUi()
  .eta <- cbind(eta.ka = c(0.1, 0.2, 0.3), eta.v = c(-0.1, -0.2, -0.3))
  .ret <- .saemEtaMatCfg(.cfg, .eta, .ui, 1:3)
  expect_equal(.ret$phiM[, 3], rep(2, 6))
  .cfg <- .saemEtaMatTestCfg()
  .ret <- .saemEtaMatCfg(.cfg, .eta[, "eta.v", drop = FALSE], .ui, 1:3, pseudoEtas = "eta.ka")
  expect_equal(.ret$phiM[, 3], rep(2, 6))
  expect_equal(.ret$phiM[, 1], 1 + .eta[rep(1:3, 2), "eta.v"])
})

test_that(".saemEtaMatCfg() errors clearly on a mismatched etaMat", {
  .cfg <- .saemEtaMatTestCfg()
  .ui <- .saemEtaMatTestUi()
  .eta <- cbind(eta.ka = c(0.1, 0.2, 0.3), eta.v = c(-0.1, -0.2, -0.3))
  expect_error(.saemEtaMatCfg(.cfg, .eta[1:2, ], .ui, 1:3), "needs 3 rows")
  expect_error(.saemEtaMatCfg(.cfg, unname(.eta[, 1, drop = FALSE]), .ui, 1:3), "needs 2 columns")
  expect_error(.saemEtaMatCfg(.cfg, cbind(.eta, eta.cl = 0), .ui, 1:3), "unknown eta.*eta.cl")
  expect_error(.saemEtaMatCfg(.cfg, .eta[, "eta.ka", drop = FALSE], .ui, 1:3), "no column for eta.*eta.v")
  .idEta <- .eta
  rownames(.idEta) <- c("1", "2", "9")
  expect_error(.saemEtaMatCfg(.cfg, .idEta, .ui, 1:3), "unknown ID.*9")
  expect_error(.saemEtaMatCfg(.cfg, .idEta[1:2, ], .ui, 1:3), "no row for ID.*3")
  expect_error(.saemEtaMatCfg(.cfg, .eta, .saemEtaMatTestUi(0.5), 1:3), "mixture")
  .old <- .uiIovEnv$iovTwoLevel
  on.exit(.uiIovEnv$iovTwoLevel <- .old)
  .uiIovEnv$iovTwoLevel <- list("eta.ka")
  expect_error(.saemEtaMatCfg(.cfg, .eta, .ui, 1:3), "IOV")
})

test_that("saemControl(etaMat=) seeds a theo_sd fit from a focei fit", {
  skip_on_cran()
  one.cmt <- function() {
    ini({
      tka <- 0.45
      tcl <- 1
      tv <- 3.45
      eta.ka ~ 0.6
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
  .f <- suppressMessages(nlmixr2(one.cmt, theo_sd, "focei", control = foceiControl(print = 0)))
  .s0 <- suppressMessages(nlmixr2(one.cmt, theo_sd, "saem", control = saemControl(print = 0)))
  .s1 <- suppressMessages(nlmixr2(one.cmt, theo_sd, "saem", control = saemControl(print = 0, etaMat = .f)))
  expect_equal(.s1$objf, .s0$objf, tolerance = 0.01)
  # IDs that are not in the data error before any iteration
  .bad <- data.frame(ID = 101:112, eta.ka = 0, eta.cl = 0, eta.v = 0)
  expect_error(
    suppressMessages(nlmixr2(one.cmt, theo_sd, "saem", control = saemControl(print = 0, etaMat = .bad))),
    "unknown ID"
  )
  expect_error(
    suppressMessages(nlmixr2(one.cmt, theo_sd, "saem", control = saemControl(print = 0, etaMat = matrix(0, 3, 3)))),
    "needs 12 rows"
  )
})
