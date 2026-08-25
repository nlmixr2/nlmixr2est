test_that("the focei-family simulation models are the default method", {
  # they used to build a `predOnly`-based model and discard it, and to call the
  # default method twice -- identical output, three times the work.  Aliases,
  # not copies, so this fails the moment a body grows back.
  for (.m in c("focei", "foce", "focep", "fo", "foi", "posthoc")) {
    expect_identical(get(paste0("getBaseSimModelFit.", .m)),
                     getBaseSimModelFit.default,
                     info = .m)
  }
})

test_that("getBaseSimModelFit() lowers the fit's model", {
  # a fit-shaped stand-in: the method only ever reads `$ui`
  .mod <- function() {
    ini({
      tka <- 0.5
      addSd <- 0.7
    })
    model({
      ka <- exp(tka)
      d/dt(depot) <- -ka * depot
      cp <- depot
      cp ~ add(addSd)
    })
  }
  .ui <- rxode2::rxode2(.mod)
  .obj <- list(ui = .ui)
  .lst <- list(.obj)
  class(.lst) <- c("focei", "getBaseSimModelFit")
  expect_equal(getBaseSimModelFit(.lst), rxode2::getBaseSimModel(.ui))
})

nmTest({
  test_that("a focei fit lowers to the same simulation model as any other fit", {
    skip_if(is.null(one.compartment.fit.focei))
    .fit <- one.compartment.fit.focei
    expect_equal(rxode2::rxNorm(.fit$simulationModel),
                 rxode2::rxNorm(eval(rxode2::getBaseSimModel(.fit$ui))))
  })
})
