test_that(".foceiPreProcessData works with data.frame and tibble data", {
  model <- function() {
    ini({
      foo <- 1
    })
    model({
      bar <- foo
      bar ~ add(foo)
    })
  }
  ui <- nlmixr(model)
  env_orig <- new.env()
  df <- data.frame(ID=1, DV=1:2, time=1:2)
  .foceiPreProcessData(data = df, env = env_orig, ui = ui)
  expect_equal(env_orig$dataSav$nlmixrRowNums, c(NA, 1, 2))
  tib <- tibble::tibble(ID=1, DV=1:2, time=1:2)
  .foceiPreProcessData(data = tib, env = env_orig, ui = ui)
  expect_equal(env_orig$dataSav$nlmixrRowNums, c(NA, 1, 2))
})
