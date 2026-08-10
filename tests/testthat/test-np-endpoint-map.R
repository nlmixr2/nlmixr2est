## Endpoint mapping for the npag/npb residual moment warm start
## (nlmixr2/nlmixr2est#856).
##
## Both halves of the mapping used to resolve "I could not match this" to endpoint
## 0, which is indistinguishable from a correct answer for the (very common)
## single-endpoint model.  A misattributed observation contributes its err^2 to the
## wrong endpoint's moment bucket, and when there is a lone scale per endpoint that
## bucket IS the estimate.  These pin the -1 ("no endpoint, drop it") sentinel.
##
## Mapping only, no fits, so this file stays in the push/PR subset.

nmTest({

  test_that("an observation cmt matching no endpoint maps to -1, not endpoint 0", {
    ## two endpoints at cmt 3 and 4 (predDf order)
    .ec <- c(3L, 4L)
    expect_equal(nlmixr2est:::npEndpointForCmt_(c(3L, 4L), .ec), c(0L, 1L))
    ## a cmt that names neither endpoint is dropped, NOT filed under the first
    expect_equal(nlmixr2est:::npEndpointForCmt_(5L, .ec), -1L)
    expect_equal(nlmixr2est:::npEndpointForCmt_(1L, .ec), -1L)
    ## getIndCmt() returns NA_INTEGER for "no compartment recorded here"
    expect_equal(nlmixr2est:::npEndpointForCmt_(NA_integer_, .ec), -1L)
    ## cmt values need not be sequential
    expect_equal(nlmixr2est:::npEndpointForCmt_(c(7L, 2L), c(2L, 7L)), c(1L, 0L))
  })

  test_that("a single-endpoint model maps every observation to endpoint 0", {
    ## a single-endpoint model carries no CMT covariate, so getIndCmt() returns 1
    ## rather than predDf$cmt (3 here) -- there is nothing to match, and every
    ## observation is that one endpoint.  This is why the fallback could not simply
    ## become -1 everywhere.
    expect_equal(nlmixr2est:::npEndpointForCmt_(c(1L, 3L, NA_integer_), 3L),
                 c(0L, 0L, 0L))
  })

  test_that("no endpoints at all drops every observation", {
    expect_equal(nlmixr2est:::npEndpointForCmt_(c(1L, 3L), integer(0)), c(-1L, -1L))
  })

  test_that("a residual parameter naming no endpoint maps to -1 and warns", {
    .endVar <- c("cp", "eff")
    expect_equal(nlmixr2est:::.npResidEndpointIdx(c("cp", "eff", "cp"), .endVar),
                 c(0L, 1L, 0L))
    expect_warning(
      .e <- nlmixr2est:::.npResidEndpointIdx(c("cp", "nosuch"), .endVar),
      "endpoint unknown")
    expect_equal(.e, c(0L, -1L))
  })

  test_that("an unavailable predDf drops every residual parameter", {
    ## the tryCatch around predDf$cond yields character(0); every match() then
    ## fails, which the old 0L coercion hid by charging them all to endpoint 0
    expect_warning(
      .e <- nlmixr2est:::.npResidEndpointIdx(c("cp", "eff"), character(0)),
      "endpoint unknown")
    expect_equal(.e, c(-1L, -1L))
  })

  test_that("no residual parameters is not a warning", {
    expect_silent(.e <- nlmixr2est:::.npResidEndpointIdx(character(0), c("cp")))
    expect_equal(.e, integer(0))
  })

})
