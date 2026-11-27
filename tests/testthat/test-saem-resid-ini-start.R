test_that("SAEM starts the residual error at its ini() values", {
  # .configsaem() placeholders, zeroed where the error model lacks the component
  .cfg <- function(resMod) {
    list(
      res.mod = resMod,
      ares = ifelse(resMod %in% c(2L, 3L, 7L, 8L), 0, 10),
      bres = ifelse(resMod %in% c(1L, 6L), 0, 1)
    )
  }
  .addProp <- function() {
    ini({
      tcl <- 1
      eta.cl ~ 0.3
      add.sd <- 0.7
      prop.sd <- 0.15
    })
    model({
      cl <- exp(tcl + eta.cl)
      ipre <- cl * time
      ipre ~ add(add.sd) + prop(prop.sd)
    })
  }
  .ui <- rxode2::rxode2(.addProp)
  .res <- .saemResIni(.cfg(.ui$saemResMod), .ui)
  expect_equal(c(.res$ares, .res$bres), c(0.7, 0.15))

  .prop <- function() {
    ini({
      tcl <- 1
      eta.cl ~ 0.3
      prop.sd <- 0.2
    })
    model({
      cl <- exp(tcl + eta.cl)
      ipre <- cl * time
      ipre ~ prop(prop.sd)
    })
  }
  .ui <- rxode2::rxode2(.prop)
  .res <- .saemResIni(.cfg(.ui$saemResMod), .ui)
  # the absent additive component stays zero
  expect_equal(c(.res$ares, .res$bres), c(0, 0.2))

  .add <- function() {
    ini({
      tcl <- 1
      eta.cl ~ 0.3
      add.sd <- 0.5
    })
    model({
      cl <- exp(tcl + eta.cl)
      ipre <- cl * time
      ipre ~ add(add.sd)
    })
  }
  .ui <- rxode2::rxode2(.add)
  .res <- .saemResIni(.cfg(.ui$saemResMod), .ui)
  expect_equal(c(.res$ares, .res$bres), c(0.5, 0))

  # a general-likelihood endpoint (res.mod 0) keeps its placeholders
  .res <- .saemResIni(list(res.mod = 0L, ares = 10, bres = 1), .ui)
  expect_equal(c(.res$ares, .res$bres), c(10, 1))
})
