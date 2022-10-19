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

test_that(".foceiPreProcessData errors with missing info", {
  env_orig <- new.env()
  df_noid <- data.frame(DV=1:2, time=1:2)
  df_nodv <- data.frame(ID=1, time=1:2)
  df_notime <- data.frame(ID=1, DV=1:2)
  expect_error(
    .foceiPreProcessData(data = df_noid, env = env_orig, ui = ui),
    regexp = "missing elements {'ID'}",
    fixed = TRUE
  )
  expect_error(
    .foceiPreProcessData(data = df_nodv, env = env_orig, ui = ui),
    regexp = "missing elements {'DV'}",
    fixed = TRUE
  )
  expect_error(
    .foceiPreProcessData(data = df_notime, env = env_orig, ui = ui),
    regexp = "missing elements {'TIME'}",
    fixed = TRUE
  )
})

test_that(".foceiPreProcessData works with data.frame and tibble data", {
  env_orig <- new.env()
  df <- data.frame(ID=1, DV=1:2, time=1:2)
  .foceiPreProcessData(data = df, env = env_orig, ui = ui)
  expect_equal(env_orig$dataSav$nlmixrRowNums, c(NA, 1, 2))
  tib <- tibble::tibble(ID=1, DV=1:2, time=1:2)
  .foceiPreProcessData(data = tib, env = env_orig, ui = ui)
  expect_equal(env_orig$dataSav$nlmixrRowNums, c(NA, 1, 2))
})
