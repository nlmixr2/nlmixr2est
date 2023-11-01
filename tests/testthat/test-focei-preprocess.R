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

test_that(".foceiPreProcessData works with ID as character or factor", {
  env_orig <- new.env()
  df <- data.frame(ID=c("A", "B"), DV=1:2, time=1:2)
  .foceiPreProcessData(data = df, env = env_orig, ui = ui)
  expect_equal(env_orig$dataSav$ID, rep(1L:2L, each = 2))
  # Confirm row order
  expect_equal(env_orig$dataSav$DV, c(NA, 1, NA, 2))

  df <- data.frame(ID=factor(c("A", "B")), DV=3:4, time=1:2)
  .foceiPreProcessData(data = df, env = env_orig, ui = ui)
  expect_equal(env_orig$dataSav$ID, rep(1L:2L, each = 2))
  # Confirm row order
  expect_equal(env_orig$dataSav$DV, c(NA, 3, NA, 4))

  df <- data.frame(ID=factor(c("B", "A")), DV=5:6, time=1:2)
  .foceiPreProcessData(data = df, env = env_orig, ui = ui)
  # IDs are assigned in the order that they are found not in factor order
  expect_equal(env_orig$dataSav$ID, rep(1L:2L, each = 2))
  # Confirm row order
  expect_equal(env_orig$dataSav$DV, c(NA, 5, NA, 6))
})
