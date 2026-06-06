## Tests that indTolRelax is stored in each control object and propagated to
## the foceiControl produced by the corresponding *ControlToFoceiControl helper.

test_that("foceiControl stores indTolRelax", {
  ctl <- foceiControl(indTolRelax = FALSE)
  expect_false(ctl$indTolRelax)
  expect_true(foceiControl()$indTolRelax)
})

test_that("saemControl stores indTolRelax and propagates to foceiControl", {
  ctl <- saemControl(indTolRelax = FALSE)
  expect_false(ctl$indTolRelax)
  env <- new.env(parent = emptyenv())
  env$saemControl <- ctl
  env$ui <- list(foceiSkipCov = NULL)
  fc <- .saemControlToFoceiControl(env, assign = FALSE)
  expect_false(fc$indTolRelax)

  env2 <- new.env(parent = emptyenv())
  env2$saemControl <- saemControl(indTolRelax = TRUE)
  env2$ui <- list(foceiSkipCov = NULL)
  fc2 <- .saemControlToFoceiControl(env2, assign = FALSE)
  expect_true(fc2$indTolRelax)
})

test_that("nlmControl stores indTolRelax and propagates to foceiControl", {
  ctl <- nlmControl(indTolRelax = FALSE)
  expect_false(ctl$indTolRelax)
  env <- new.env(parent = emptyenv())
  env$nlmControl <- ctl
  env$ui <- list(foceiSkipCov = NULL)
  fc <- .nlmControlToFoceiControl(env, assign = FALSE)
  expect_false(fc$indTolRelax)
})

test_that("bobyqaControl stores indTolRelax and propagates to foceiControl", {
  ctl <- bobyqaControl(indTolRelax = FALSE)
  expect_false(ctl$indTolRelax)
  env <- new.env(parent = emptyenv())
  env$bobyqaControl <- ctl
  env$ui <- list(foceiSkipCov = NULL)
  fc <- .bobyqaControlToFoceiControl(env, assign = FALSE)
  expect_false(fc$indTolRelax)
})

test_that("nlminbControl stores indTolRelax and propagates to foceiControl", {
  ctl <- nlminbControl(indTolRelax = FALSE)
  expect_false(ctl$indTolRelax)
  env <- new.env(parent = emptyenv())
  env$nlminbControl <- ctl
  env$ui <- list(foceiSkipCov = NULL)
  fc <- .nlminbControlToFoceiControl(env, assign = FALSE)
  expect_false(fc$indTolRelax)
})

test_that("n1qn1Control stores indTolRelax and propagates to foceiControl", {
  ctl <- n1qn1Control(indTolRelax = FALSE)
  expect_false(ctl$indTolRelax)
  env <- new.env(parent = emptyenv())
  env$n1qn1Control <- ctl
  env$ui <- list(foceiSkipCov = NULL)
  fc <- .n1qn1ControlToFoceiControl(env, assign = FALSE)
  expect_false(fc$indTolRelax)
})

test_that("lbfgsb3cControl stores indTolRelax and propagates to foceiControl", {
  ctl <- lbfgsb3cControl(indTolRelax = FALSE)
  expect_false(ctl$indTolRelax)
  env <- new.env(parent = emptyenv())
  env$lbfgsb3cControl <- ctl
  env$ui <- list(foceiSkipCov = NULL)
  fc <- .lbfgsb3cControlToFoceiControl(env, assign = FALSE)
  expect_false(fc$indTolRelax)
})

test_that("optimControl stores indTolRelax and propagates to foceiControl", {
  ctl <- optimControl(indTolRelax = FALSE)
  expect_false(ctl$indTolRelax)
  env <- new.env(parent = emptyenv())
  env$optimControl <- ctl
  env$ui <- list(foceiSkipCov = NULL)
  fc <- .optimControlToFoceiControl(env, assign = FALSE)
  expect_false(fc$indTolRelax)
})

test_that("nlsControl stores indTolRelax and propagates to foceiControl", {
  ctl <- nlsControl(indTolRelax = FALSE)
  expect_false(ctl$indTolRelax)
  env <- new.env(parent = emptyenv())
  env$nlsControl <- ctl
  env$ui <- list(foceiSkipCov = NULL)
  fc <- .nlsControlToFoceiControl(env, assign = FALSE)
  expect_false(fc$indTolRelax)
})

test_that("uobyqaControl stores indTolRelax and propagates to foceiControl", {
  ctl <- uobyqaControl(indTolRelax = FALSE)
  expect_false(ctl$indTolRelax)
  env <- new.env(parent = emptyenv())
  env$uobyqaControl <- ctl
  env$ui <- list(foceiSkipCov = NULL)
  fc <- .uobyqaControlToFoceiControl(env, assign = FALSE)
  expect_false(fc$indTolRelax)
})

test_that("newuoaControl stores indTolRelax and propagates to foceiControl", {
  ctl <- newuoaControl(indTolRelax = FALSE)
  expect_false(ctl$indTolRelax)
  env <- new.env(parent = emptyenv())
  env$newuoaControl <- ctl
  env$ui <- list(foceiSkipCov = NULL)
  fc <- .newuoaControlToFoceiControl(env, assign = FALSE)
  expect_false(fc$indTolRelax)
})

test_that("nlmeControlToFoceiControl always passes indTolRelax=TRUE", {
  env <- new.env(parent = emptyenv())
  env$nlmeControl <- nlmeControl()
  env$etaMat <- NULL
  env$ui <- list(foceiSkipCov = NULL)
  fc <- .nlmeControlToFoceiControl(env, assign = FALSE)
  expect_true(fc$indTolRelax)
})

test_that("foceControl inherits indTolRelax from foceiControl default", {
  ctl <- foceControl()
  expect_true(ctl$indTolRelax)
  ctl2 <- foceControl(indTolRelax = FALSE)
  expect_false(ctl2$indTolRelax)
})

test_that("laplaceControl inherits indTolRelax from foceiControl default", {
  ctl <- laplaceControl()
  expect_true(ctl$indTolRelax)
  ctl2 <- laplaceControl(indTolRelax = FALSE)
  expect_false(ctl2$indTolRelax)
})

test_that("agqControl inherits indTolRelax from foceiControl default", {
  ctl <- agqControl()
  expect_true(ctl$indTolRelax)
  ctl2 <- agqControl(indTolRelax = FALSE)
  expect_false(ctl2$indTolRelax)
})
