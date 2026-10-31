test_that('compiled FOCEi cache keys distinguish rxode2 versions and builds', {
  model <- function() {
    ini({ level <- .2; error <- .2; eta ~ .3 })
    model({ prediction <- level + eta; prediction ~ add(error) })
  }
  ui <- list(rxode2::.copyUi(suppressMessages(nlmixr2(model))))
  original <- rxUiGet.foceiModelDigest(ui)
  expect_identical(rxUiGet.foceiModelDigest(ui), original)
  for (field in c('version', 'md5')) {
    other <- rxode2::rxVersion()
    other[field] <- paste0(other[field], '-other')
    updated <- with_mocked_bindings(
      rxUiGet.foceiModelDigest(ui),
      rxVersion = function(...) other,
      .package = 'rxode2'
    )
    expect_false(identical(updated, original), info = field)
  }
})
