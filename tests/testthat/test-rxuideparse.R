test_that("foceiControl() deparse", {

  expect_equal(rxUiDeparse.foceiControl(foceiControl(innerOpt="BFGS", scaleType="norm", normType="std", derivMethod="central", covDerivMethod="forward", covMethod="s",diagXform="identity", addProp= "combined1", eventType="forward", optimHessType="forward"), "ctl"),
               quote(ctl <- foceiControl(derivMethod = "central", covDerivMethod = "forward",
                                         covMethod = "s", diagXform = "identity", optimHessType = "forward",
                                         innerOpt = "BFGS", scaleType = "norm", normType = "std",
                                         eventType = "forward", addProp = "combined1")))

  expect_equal(rxUiDeparse.foceiControl(foceiControl(eventType="forward"), "ctl"),
               quote(ctl <- foceiControl(eventType = "forward")))

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

test_that("bobyqaControl()",{

  expect_equal(rxUiDeparse.bobyqaControl(bobyqaControl(), "var"),
               quote(var <- bobyqaControl()))

  expect_equal(rxUiDeparse.bobyqaControl(bobyqaControl(scaleType="multAdd"), "var"),
               quote(var <- bobyqaControl(scaleType = "multAdd")))

})

test_that("lbfgsb3cControl()", {
  expect_equal(rxUiDeparse.lbfgsb3cControl(lbfgsb3cControl(), "var"),
               quote(var <- lbfgsb3cControl()))

  expect_equal(rxUiDeparse.lbfgsb3cControl(lbfgsb3cControl(normType="len"), "var"),
               quote(var <- lbfgsb3cControl(normType = "len")))
})

test_that("n1qn1Control()", {
  expect_equal(rxUiDeparse.n1qn1Control(n1qn1Control(), "var"),
               quote(var <- n1qn1Control()))

  expect_equal(rxUiDeparse.n1qn1Control(n1qn1Control(covMethod="n1qn1"), "var"),
               quote(var <- n1qn1Control(covMethod = "n1qn1")))
})

test_that("newuoaControl()", {
  expect_equal(rxUiDeparse.newuoaControl(newuoaControl(), "var"),
               quote(var <- newuoaControl()))
  expect_equal(rxUiDeparse.newuoaControl(newuoaControl(addProp="combined1"), "var"),
               quote(var <- newuoaControl(addProp = "combined1")))
})

test_that("nlmeControl()", {
  expect_equal(rxUiDeparse.nlmeControl(nlmeControl(), "var"),
               quote(var <- nlmeControl()))

  expect_equal(rxUiDeparse.nlmeControl(nlmeControl(opt="nlm"), "var"),
               quote(var <- nlmeControl(opt = "nlm")))
})

test_that("nlminbControl()", {
  expect_equal(rxUiDeparse.nlminbControl(nlminbControl(), "var"),
               quote(var <- nlminbControl()))

  expect_equal(rxUiDeparse.nlminbControl(nlminbControl(solveType="grad"), "var"),
               quote(var <- nlminbControl(solveType = "grad")))
})

test_that("nlmControl()", {
  expect_equal(rxUiDeparse.nlmControl(nlmControl(), "var"),
               quote(var <- nlmControl()))
  expect_equal(rxUiDeparse.nlmControl(nlmControl(covMethod="r"), "var"),
               quote(var <- nlmControl(covMethod = "r")))
})

test_that("nlsControl()", {
  expect_equal(rxUiDeparse.nlsControl(nlsControl(), "var"),
               quote(var <- nlsControl()))
  expect_equal(rxUiDeparse.nlsControl(nlsControl(algorithm="port"), "var"),
               quote(var <- nlsControl(algorithm = "port")))
})

test_that("optimControl()", {
  expect_equal(rxUiDeparse.optimControl(optimControl(), "var"),
               quote(var <- optimControl()))
  expect_equal(rxUiDeparse.optimControl(optimControl(method="L-BFGS-B", covMethod = "optim"), "var"),
               quote(var <- optimControl(method = "L-BFGS-B", covMethod="optim")))

  expect_equal(rxUiDeparse.optimControl(optimControl(eventType="forward"), "var"),
               quote(var <- optimControl(eventType = "forward")))
})

test_that("uobyqaControl()", {
  expect_equal(rxUiDeparse.uobyqaControl(uobyqaControl(), "var"),
               quote(var <- uobyqaControl()))
  expect_equal(rxUiDeparse.uobyqaControl(uobyqaControl(scaleTo=4), "var"),
               quote(var <- uobyqaControl(scaleTo = 4)))
})


test_that("tableControl()", {
  expect_equal(rxUiDeparse.tableControl(tableControl(), "var"),
               quote(var <- tableControl()))
  expect_equal(rxUiDeparse.tableControl(tableControl(censMethod="epred"), "var"),
               quote(var <- tableControl(censMethod = "epred")))
})
