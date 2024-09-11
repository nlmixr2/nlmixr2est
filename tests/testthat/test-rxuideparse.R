## S3method(rxUiDeparse,foceiControl)
test_that("foceiControl() deparse", {

  expect_equal(rxUiDeparse(foceiControl(innerOpt="BFGS", scaleType="norm", normType="std", derivMethod="central", covDerivMethod="forward", covMethod="s",diagXform="identity", addProp= "combined1", rxControl=rxControl(returnType="data.frame")), "ctl"),

               quote(ctl <- foceiControl(derivMethod = "central", covDerivMethod = "forward",
                                         covMethod = "s", diagXform = "identity", innerOpt = "BFGS",
                                         scaleType = "norm", normType = "std", addProp = "combined1",
                                         rxControl = rxControl(returnType = "data.frame"))))

})

## S3method(rxUiDeparse,n1qn1Control)
test_that("n1qn1Control()", {
  expect_equal(rxUiDeparse(n1qn1Control(), "ctl"),
               quote(ctl <- n1qn1Control()))
  expect_equal(rxUiDeparse(n1qn1Control(nsim=100), "ctl"),
               quote(ctl <- n1qn1Control(nsim=100)))
})

## S3method(rxUiDeparse,nlminbControl)
test_that("nlminbControl()", {
  expect_equal(rxUiDeparse(nlminbControl(), "ctl"),
               quote(ctl <- nlminbControl()))
  expect_equal(rxUiDeparse(nlminbControl(eval.max=100), "ctl"),
               quote(ctl <- nlminbControl(eval.max=100)))
})

## S3method(rxUiDeparse,newuoaControl)
test_that("newuoaControl()", {
  expect_equal(rxUiDeparse(newuoaControl(), "ctl"),
               quote(ctl <- newuoaControl()))
  expect_equal(rxUiDeparse(newuoaControl(maxfun=100), "ctl"),
               quote(ctl <- newuoaControl(maxfun=100)))
})
## S3method(rxUiDeparse,nlmControl)
test_that("nlmControl()", {
  expect_equal(rxUiDeparse(nlmControl(), "ctl"),
               quote(ctl <- nlmControl()))
  expect_equal(rxUiDeparse(nlmControl(fscale=2), "ctl"),
               quote(ctl <- nlmControl(fscale=2)))
})
## S3method(rxUiDeparse,nlmeControl)
test_that("nlmeControl()", {
  expect_equal(rxUiDeparse(nlmeControl(), "ctl"),
               quote(ctl <- nlmeControl()))
  expect_equal(rxUiDeparse(nlmeControl(maxIter=1000), "ctl"),
               quote(ctl <- nlmeControl(maxIter=1000)))
})
## S3method(rxUiDeparse,nlsControl)
test_that("nlsControl()", {
  expect_equal(rxUiDeparse(nlsControl(), "ctl"),
               quote(ctl <- nlsControl()))
  expect_equal(rxUiDeparse(nlsControl(maxiter=2), "ctl"),
               quote(ctl <- nlsControl(maxiter=2)))
})
## S3method(rxUiDeparse,optimControl)
test_that("optimControl()", {
  expect_equal(rxUiDeparse(optimControl(), "ctl"),
               quote(ctl <- optimControl()))
  expect_equal(rxUiDeparse(optimControl(method="L-BFGS-B"), "ctl"),
               quote(ctl <- optimControl(method = "L-BFGS-B", covMethod="optim")))
})

## S3method(rxUiDeparse,saemControl)
test_that("saemControl()", {
  expect_equal(rxUiDeparse(saemControl(), "ctl"),
               quote(ctl <- saemControl()))
  expect_equal(rxUiDeparse(saemControl(nEm=4), "ctl"),
               quote(ctl <- saemControl(nEm=4)))
})

## S3method(rxUiDeparse,tableControl)
test_that("tableControl()", {
  expect_equal(rxUiDeparse(tableControl(), "ctl"),
               quote(ctl <- tableControl()))
  expect_equal(rxUiDeparse(tableControl(censMethod="epred", npde=TRUE), "ctl"),
               quote(ctl <- tableControl(npde = TRUE, censMethod = "epred")))
})
