test_that("foceiControl() deparse", {

  expect_equal(rxUiDeparse(foceiControl(innerOpt="BFGS", scaleType="norm", normType="std", derivMethod="central", covDerivMethod="forward", covMethod="s",diagXform="identity", addProp= "combined1", rxControl=rxControl(returnType="data.frame")), "ctl"),

               quote(ctl <- foceiControl(derivMethod = "central", covDerivMethod = "forward",
                                         covMethod = "s", diagXform = "identity", innerOpt = "BFGS",
                                         scaleType = "norm", normType = "std", addProp = "combined1",
                                         rxControl = rxControl(returnType = "data.frame"))))

})

test_that("n1qn1Control()", {
  expect_equal(rxUiDeparse(n1qn1Control(), "ctl"),
               quote(ctl <- n1qn1Control()))
  expect_equal(rxUiDeparse(n1qn1Control(nsim=100), "ctl"),
               quote(ctl <- n1qn1Control(nsim=100)))
})

test_that("nlminbControl()", {
  expect_equal(rxUiDeparse(nlminbControl(), "ctl"),
               quote(ctl <- nlminbControl()))
  expect_equal(rxUiDeparse(nlminbControl(eval.max=100), "ctl"),
               quote(ctl <- nlminbControl(eval.max=100)))
})
