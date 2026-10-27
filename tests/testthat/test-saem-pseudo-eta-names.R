# saem's temporary etas are found by its record of them (the transform specs
# with pseudoEta = TRUE), not by the rx.eta.<name> pattern, which an occasion
# eta eta.<x> also takes in the IOV rewrite (rx.eta.<x>.<occ>)

# a stub fit whose environment has a ranef, a shrink and `transforms` as the
# temporary etas' record
.pseudoEtaStub <- function(transforms) {
  .env <- new.env(parent = emptyenv())
  .env$ranef <- data.frame(ID = 1:2, eta.cl = c(0.1, 0.2), rx.eta.ka.occ.1 = c(0.3, 0.4), rx.eta.lam = c(0.5, 0.6))
  .env$shrink <- data.frame(eta.cl = 1, rx.eta.ka.occ.1 = 2, rx.eta.lam = 3, row.names = "sd")
  .env$saemPseudoTransforms <- transforms
  list(env = .env)
}

test_that(".saemPseudoTransformNames() gives the thetas of the specs marked pseudoEta", {
  .tr <- list(
    list(name = "lam", internalName = "rxBoundedTr.lam", pseudoEta = TRUE),
    list(name = "add.sd", internalName = "rxBoundedTr.add.sd"),
    list(name = "nu", internalName = "rxBoundedTr.nu", pseudoEta = FALSE)
  )
  expect_identical(.saemPseudoTransformNames(.tr), "lam")
  expect_identical(.saemPseudoTransformNames(NULL), character(0))
  expect_identical(.saemPseudoTransformNames(list()), character(0))
})

test_that("the temporary-eta finalize drops only the recorded etas", {
  # no record (every estimator but saem with temporary etas): an occasion eta
  # named rx.eta.<x>.<occ> is kept
  .r <- .saemPseudoEtaFinalize(.pseudoEtaStub(NULL))
  expect_identical(names(.r$env$ranef), c("ID", "eta.cl", "rx.eta.ka.occ.1", "rx.eta.lam"))
  expect_identical(names(.r$env$shrink), c("eta.cl", "rx.eta.ka.occ.1", "rx.eta.lam"))
  # a record for lam drops rx.eta.lam alone
  .r <- .saemPseudoEtaFinalize(.pseudoEtaStub(list(list(
    name = "lam",
    internalName = "rxBoundedTr.lam",
    pseudoEta = TRUE
  ))))
  expect_identical(names(.r$env$ranef), c("ID", "eta.cl", "rx.eta.ka.occ.1"))
  expect_identical(names(.r$env$shrink), c("eta.cl", "rx.eta.ka.occ.1"))
})
