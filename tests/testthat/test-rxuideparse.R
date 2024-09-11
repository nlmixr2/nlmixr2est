test_that("foceiControl() deparse", {

  expect_equal(rxUiDeparse.foceiControl(foceiControl(innerOpt="BFGS", scaleType="norm", normType="std", derivMethod="central", covDerivMethod="forward", covMethod="s",diagXform="identity", addProp= "combined1"), "ctl"),
               quote(ctl <- foceiControl(derivMethod = "central", covDerivMethod = "forward",
                                         covMethod = "s", diagXform = "identity", innerOpt = "BFGS",
                                         scaleType = "norm", normType = "std", addProp = "combined1")))

})
