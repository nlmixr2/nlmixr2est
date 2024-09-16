test_that("foceiControl() deparse", {

  expect_equal(rxUiDeparse.foceiControl(foceiControl(innerOpt="BFGS", scaleType="norm", normType="std", derivMethod="central", covDerivMethod="forward", covMethod="s",diagXform="identity", addProp= "combined1"), "ctl"),
               quote(ctl <- foceiControl(derivMethod = "central", covDerivMethod = "forward",
                                         covMethod = "s", diagXform = "identity", innerOpt = "BFGS",
                                         scaleType = "norm", normType = "std", addProp = "combined1")))

  expect_equal(rxUiDeparse.foceiControl(foceiControl(), "ctl"),
               quote(ctl <- foceiControl()))

  expect_warning(rxUiDeparse.foceiControl(foceiControl(outerOpt=optim), "ctl"),
                 "reset")

})

test_that("saemControl() deparse", {

  expect_equal(rxUiDeparse.saemControl(saemControl(nBurn=2, nEm=2, nmc=7, nu=c(3, 3, 3)),
                                       "ctl"),
               quote(ctl <- saemControl(nBurn = 2, nEm = 2, nmc = 7,
                                        nu = c(3, 3, 3))))
  expect_equal(rxUiDeparse.saemControl(saemControl(), "ctl"),
               quote(ctl <- saemControl()))
})
