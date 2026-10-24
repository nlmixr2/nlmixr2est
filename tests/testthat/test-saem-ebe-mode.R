# saemControl(ebe="mode") replaces SAEM's approximate conditional-mean etas with
# each subject's conditional mode (the Monolix "EBEs" task), computed in C++
# after the last SAEM iteration.  The SAEM etas stay available as fit$etaSaem.

test_that("saem ebe='mode' returns the conditional mode", {
  skip_on_cran()

  # nolint start: object_usage_linter.
  .mod <- function() {
    ini({
      tka <- 0.45
      tcl <- 1
      tv <- 3.45
      add.sd <- 0.7
      eta.ka ~ 0.6
      eta.cl ~ 0.3
      eta.v ~ 0.1
    })
    model({
      ka <- exp(tka + eta.ka)
      cl <- exp(tcl + eta.cl)
      v <- exp(tv + eta.v)
      linCmt() ~ add(add.sd)
    })
  }
  # nolint end

  .ctl <- saemControl(print = 0, nBurn = 50, nEm = 50)
  fit <- .nlmixr(.mod, nlmixr2data::theo_sd, est = "saem", control = .ctl)
  expect_equal(fit$ebe, "mode")
  .es <- fit$etaSaem
  expect_true(is.data.frame(.es))
  expect_equal(names(.es), c("ID", "eta.ka", "eta.cl", "eta.v"))
  .etaCols <- c("eta.ka", "eta.cl", "eta.v")
  .mode <- as.matrix(fit$eta[, .etaCols])
  .saem <- as.matrix(.es[, .etaCols])
  expect_false(isTRUE(all.equal(.mode, .saem, check.attributes = FALSE)))

  # FOCEi's inner problem at the same final parameters is the same conditional
  # density, so its converged etas are the mode too
  .f2 <- .nlmixr(fit, est = "focei",
                 control = foceiControl(print = 0, maxOuterIterations = 0L,
                                        covMethod = ""))
  expect_equal(as.matrix(.f2$eta[, .etaCols]), .mode,
               tolerance = 1e-3, ignore_attr = TRUE)

  .fitS <- .nlmixr(.mod, nlmixr2data::theo_sd, est = "saem",
                   control = saemControl(print = 0, nBurn = 50, nEm = 50,
                                         ebe = "saem"))
  expect_equal(.fitS$ebe, "saem")
  expect_equal(as.matrix(.fitS$eta[, .etaCols]),
               as.matrix(.fitS$etaSaem[, .etaCols]), ignore_attr = TRUE)
  # the SAEM run itself does not depend on the ebe choice
  expect_equal(.fitS$etaSaem, fit$etaSaem)
  expect_equal(.fitS$theta, fit$theta)

  # each subject's search is independent, so threads do not change the answer
  .fit2 <- .nlmixr(.mod, nlmixr2data::theo_sd, est = "saem",
                   control = saemControl(print = 0, nBurn = 50, nEm = 50,
                                         rxControl = rxode2::rxControl(cores = 2L)))
  expect_equal(as.matrix(.fit2$eta[, .etaCols]), .mode, ignore_attr = TRUE)
})

test_that("saem ebe='mode' handles a single eta (1-D search)", {
  skip_on_cran()

  # nolint start: object_usage_linter.
  .mod <- function() {
    ini({
      tka <- 0.45
      tcl <- 1
      tv <- 3.45
      add.sd <- 0.7
      eta.cl ~ 0.3
    })
    model({
      ka <- exp(tka)
      cl <- exp(tcl + eta.cl)
      v <- exp(tv)
      linCmt() ~ add(add.sd)
    })
  }
  # nolint end

  fit <- .nlmixr(.mod, nlmixr2data::theo_sd, est = "saem",
                 control = saemControl(print = 0, nBurn = 50, nEm = 50))
  .f2 <- .nlmixr(fit, est = "focei",
                 control = foceiControl(print = 0, maxOuterIterations = 0L,
                                        covMethod = ""))
  expect_equal(.f2$eta$eta.cl, fit$eta$eta.cl, tolerance = 1e-3)
  expect_false(isTRUE(all.equal(fit$eta$eta.cl, fit$etaSaem$eta.cl)))
})
