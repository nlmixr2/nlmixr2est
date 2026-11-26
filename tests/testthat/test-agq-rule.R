test_that("Gauss-Legendre nodes match the closed form", {
  .gl <- .gaussLegendre(3L)
  expect_equal(.gl$x, c(-sqrt(3 / 5), 0, sqrt(3 / 5)))
  expect_equal(.gl$w, c(5, 8, 5) / 9)
})

test_that("normal-measure rules integrate constants and second moments", {
  .moments <- function(g) {
    z <- g$x * sqrt(2)
    p <- exp(rowSums(log(g$w)))
    list(sum = sum(p), m2 = crossprod(z * sqrt(p)), m1 = colSums(z * p))
  }
  for (neta in 1:4) {
    .s <- .moments(.agqStieltjes(neta))
    expect_equal(.s$sum, 1)
    expect_equal(.s$m1, rep(0, neta))
    expect_equal(.s$m2, diag(neta), ignore_attr = TRUE)
    .l <- .moments(.agqLegendre(neta, 3L))
    expect_equal(.l$sum, 1)
    expect_equal(.l$m1, rep(0, neta))
  }
})

test_that("Stieltjes grid size and tail rule", {
  # 5 regions (centre, 3 shells, tail) x 2*ns*neta*(neta-1) directions
  expect_equal(.agqStieltjes(3)$n, 5 * 24)
  expect_equal(.agqStieltjes(1)$n, 5 * 2)
  expect_equal(.agqStieltjes(2, c(3, 1, 0.6, 0.9))$n, 5 * 4)
  # r1 > 0.9999 drops the tail region
  expect_equal(.agqStieltjes(2, c(3, 1, 0.6, 0.99999))$n, 4 * 4)
  expect_false(.agqStieltjes(2)$first)
  expect_true(.agqLegendre(2, 3L)$first)
  expect_equal(.agqLegendre(2, 3L)$x[1, ], c(0, 0))
})

test_that("agqControl validates the rule options", {
  expect_equal(agqControl()$agqRule, "hermite")
  expect_equal(agqControl(agqRule = "stieltjes")$agqGrid, c(nr = 3, ns = 2, r0 = 0.6, r1 = 0.9))
  expect_error(agqControl(agqGrid = c(3, 1, 0.9, 0.6)), "agqGrid")
  expect_error(agqControl(agqGrid = c(0, 1, 0.6, 0.9)), "agqGrid")
  expect_error(agqControl(agqNsd = 0), "agqNsd")
  expect_error(foceiControl(agqRule = "legendre"), "nAGQ > 1")
  expect_error(.agqStieltjes(3, maxNodes = 10), "exceeds the limit")
})

test_that("every rule is exact for a model linear in eta", {
  skip_on_cran()
  set.seed(42)
  d <- do.call(
    rbind,
    lapply(1:20, function(i) {
      e <- rnorm(2, 0, c(0.7, 0.5))
      t <- c(1, 2, 4, 8)
      data.frame(id = i, time = t, dv = 1 + e[1] + (0.5 + e[2]) * t + rnorm(4, 0, 0.3))
    })
  )
  f <- function() {
    ini({
      ta <- 1
      tb <- 0.5
      eta.a ~ 0.5
      eta.b ~ 0.2
      add.sd <- 0.3
    })
    model({
      a <- ta + eta.a
      b <- tb + eta.b
      y <- a + b * time
      y ~ add(add.sd)
    })
  }
  .obj <- function(...) {
    .fit <- suppressMessages(nlmixr2(
      f,
      d,
      est = "agq",
      control = agqControl(maxOuterIterations = 0, covMethod = "", print = 0, ...)
    ))
    list(obj = .fit$objective, row = rownames(.fit$objDf))
  }
  .lap <- .obj(nAGQ = 1)$obj
  .st <- .obj(agqRule = "stieltjes")
  expect_equal(.st$obj, .lap, tolerance = 1e-6)
  expect_equal(.st$row, "AGQ-Stieltjes")
  .gl <- .obj(nAGQ = 4, agqRule = "legendre")
  expect_equal(.gl$obj, .lap, tolerance = 1e-6)
  expect_equal(.gl$row, "AGQ-GL4")
})

test_that("the rules converge to the same nonlinear integral", {
  skip_on_cran()
  p <- pump
  p$dv <- p$y
  p$time <- 0
  f <- function() {
    ini({
      t1 <- 1
      t2 <- 1
      t3 <- 1
      t4 <- 1
      eta1 ~ 1
    })
    model({
      if (group == 1) {
        lp <- t1 + t2 * logtstd
      } else {
        lp <- t3 + t4 * logtstd
      }
      lp <- lp + eta1
      lam <- exp(lp)
      y ~ dpois(lam)
    })
  }
  .obj <- function(...) {
    suppressMessages(nlmixr2(
      f,
      p,
      est = "agq",
      control = agqControl(maxOuterIterations = 0, covMethod = "", print = 0, ...)
    ))$objective
  }
  .ref <- .obj(nAGQ = 9)
  .st <- .obj(agqRule = "stieltjes", agqGrid = c(30, 1, 0.2, 0.99999))
  .gl <- .obj(nAGQ = 15, agqRule = "legendre", agqNsd = 5)
  expect_lt(abs(.st - .ref), 0.05)
  expect_lt(abs(.gl - .ref), 0.01)
  # an estimation run reaches the Gauss-Hermite estimates
  .fit <- suppressMessages(suppressWarnings(nlmixr2(
    f,
    p,
    est = "agq",
    control = agqControl(agqRule = "stieltjes", print = 0)
  )))
  expect_equal(unname(fixef(.fit)[1:4]), c(2.96, -0.426, 1.80, 0.611), tolerance = 0.02)
})
