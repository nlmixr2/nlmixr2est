
test_that("test timing engine", {

  env <- new.env(parent=emptyenv())
  env$time <- data.frame(t0=1)

  nlmixrWithTiming("time1", {
    Sys.sleep(0.25)
    # note this can be nested, time1 will exclude the timing from time2
    nlmixrWithTiming("time2", {
      Sys.sleep(0.25)
      nlmixrWithTiming("t0", {
        Sys.sleep(0.1)
      }, envir=env)
    }, envir=env)
  }, envir=env)

  expect_true(env$time$time2 < 0.4)
  expect_true(env$time$time1 < 0.4)
  expect_true(env$time$t0 > 1)

  nlmixrAddTiming(env, "time2", 3)

  expect_true(env$time$time2 > 3)

  nlmixrAddTiming(env, "new", 3.5)

  expect_true(env$time$new > 3)

})
