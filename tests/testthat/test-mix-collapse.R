test_that(".mixCollapseCheck is silent for well-separated classes", {
  post <- cbind(c(0.99, 0.02, 0.9), c(0.01, 0.98, 0.1))
  ll <- cbind(c(-10, -30, -12), c(-15, -12, -14))
  expect_silent(r <- .mixCollapseCheck(post, ll))
  expect_length(r$empty, 0L)
  expect_equal(nrow(r$same), 0L)
})

test_that(".mixCollapseCheck names an empty class", {
  post <- cbind(c(0.995, 0.999, 0.992), c(0.005, 0.001, 0.008))
  expect_warning(
    r <- .mixCollapseCheck(post),
    "mixture class 2 has no subjects \\(max posterior < 0.01\\); try new starts"
  )
  expect_equal(r$empty, 2L)
  # the boundary itself is not empty
  post[1, ] <- c(0.99, 0.01)
  expect_silent(.mixCollapseCheck(post))
  expect_warning(.mixCollapseCheck(post, emptyTol = 0.0100001), "class 2")
})

test_that(".mixCollapseCheck names an indistinguishable pair", {
  post <- cbind(c(0.9, 0.05, 0.2), c(0.06, 0.03, 0.5), c(0.04, 0.92, 0.3))
  ll <- cbind(c(-10, -30, -12), c(-30, -12, -12.5), c(-30.01, -12.02, -12.49))
  expect_warning(
    r <- .mixCollapseCheck(post, ll),
    "mixture classes 2 and 3 are indistinguishable at the estimates"
  )
  expect_equal(r$same, matrix(c(2L, 3L), ncol = 2L))
  # one subject that tells them apart is enough
  ll[3, 3] <- -12.5 + 0.05
  expect_silent(.mixCollapseCheck(post, ll))
  ll[3, 3] <- -12.5 + 0.0499
  expect_warning(.mixCollapseCheck(post, ll), "classes 2 and 3")
})

test_that(".mixCollapseCheck reports an empty class once, not also as a pair", {
  post <- cbind(c(0.999, 0.998), c(0.001, 0.002))
  ll <- cbind(c(-10, -11), c(-10, -11))
  w <- character(0)
  withCallingHandlers(
    .mixCollapseCheck(post, ll),
    warning = function(e) {
      w <<- c(w, conditionMessage(e))
      invokeRestart("muffleWarning")
    }
  )
  expect_equal(w, "mixture class 2 has no subjects (max posterior < 0.01); try new starts")
})

test_that(".mixCollapseCheck ignores subjects with no evidence for either class", {
  post <- cbind(c(1, 0, 0), c(0, 0.5, 0.5), c(0, 0.5, 0.5))
  # subject 1 belongs to class 1: both other log-posteriors are -Inf
  ll <- log(post) - rep(log(c(1 / 3, 1 / 3, 1 / 3)), each = 3)
  expect_warning(.mixCollapseCheck(post, ll), "classes 2 and 3")
  expect_silent(.mixCollapseCheck(post, NULL))
  expect_silent(.mixCollapseCheck(post[, 1, drop = FALSE]))
})

test_that("mixture warnings stay under 75 characters", {
  expect_lt(nchar("mixture class 10 has no subjects (max posterior < 0.01); try new starts"), 75)
  expect_lt(nchar("mixture classes 10 and 11 are indistinguishable at the estimates"), 75)
})

nmTest({
  test_that("a focei mixture fit warns on an empty or merged class, not on separated ones", {
    rxode2::rxWithSeed(42, {
      clSim <- ifelse(rbinom(20, 1, 0.5) == 1, 1.2, 6.0)
      times <- c(0.5, 1, 2, 4, 8, 12, 24)
      d <- do.call(
        rbind,
        lapply(seq_along(clSim), function(i) {
          k <- clSim[i] / 24
          cp <- 100 * 1.5 / (24 * (1.5 - k)) * (exp(-k * times) - exp(-1.5 * times)) + rnorm(length(times), 0, 0.05)
          data.frame(
            ID = i,
            TIME = c(0, times),
            AMT = c(100, rep(0, 7)),
            EVID = c(1, rep(0, 7)),
            DV = c(0, pmax(cp, 0)),
            CMT = c(1, rep(2, 7))
          )
        })
      )
    })
    mixMod <- function() {
      ini({
        tka <- log(1.5)
        tcl1 <- fix(log(1.2))
        tcl2 <- fix(log(6))
        tv <- log(24)
        p1 <- 0.5
        eta.cl ~ 0.01
        eta.v ~ 0.01
        add.sd <- 0.05
      })
      model({
        ka <- exp(tka)
        cl <- mix(exp(tcl1 + eta.cl), p1, exp(tcl2 + eta.cl))
        v <- exp(tv + eta.v)
        d / dt(depot) <- -ka * depot
        d / dt(center) <- ka * depot - cl / v * center
        cp <- center / v
        cp ~ add(add.sd)
      })
    }
    .ctl <- foceiControl(print = 0, maxOuterIterations = 0L, covMethod = "", calcTables = FALSE)
    .mixWarn <- function(fit) grep("^mixture class", fit$runInfo, value = TRUE)

    fit <- suppressMessages(nlmixr2(mixMod, d, "focei", control = .ctl))
    expect_length(.mixWarn(fit), 0L)

    fit <- suppressMessages(nlmixr2(ini(mixMod, tcl2 = log(1.2)), d, "focei", control = .ctl))
    expect_equal(.mixWarn(fit), "mixture classes 1 and 2 are indistinguishable at the estimates")

    fit <- suppressMessages(nlmixr2(ini(mixMod, tcl2 = log(500)), d, "focei", control = .ctl))
    expect_equal(.mixWarn(fit), "mixture class 2 has no subjects (max posterior < 0.01); try new starts")
  })
})
