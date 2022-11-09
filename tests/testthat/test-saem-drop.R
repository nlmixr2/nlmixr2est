nmTest({
  test_that("saem dropping parameters", {
    .nm <- loadNamespace("nlmixr2est")

    muRefDataFrame <-
      structure(list(theta = c("tka", "tcl", "tv"),
                     eta = c("eta.ka", "eta.cl", "eta.v"),
                     level = c("id", "id", "id")),
                row.names = c(NA, -3L),
                class = "data.frame")

    muRefCovariateDataFrame <-
      structure(list(theta = c("tcl", "tcl", "tv", "tv", "tv", "tvp", "tvp", "tvp"),
                     covariate = c("age", "sex", "age", "sex", "wt", "age", "sex", "wt"),
                     covariateParameter = c("cl.age", "cl.sex", "v.age", "v.sex", "v.wt", "vp.age", "vp.sex", "vp.wt")),
                row.names = c(NA, -8L),
                class = "data.frame")

    expect_equal(.nm$.saemDropParameters(quote(ka <- exp(tka + eta.ka)), muRefDataFrame, muRefCovariateDataFrame),
                 quote(ka <- exp(tka)))

    expect_equal(.nm$.saemDropParameters(quote(ka <- exp(eta.ka + tka)), muRefDataFrame, muRefCovariateDataFrame),
                 quote(ka <- exp(tka)))


    expect_equal(.nm$.saemDropParameters(quote(ka <- exp(eta.ka + tka + 1)), muRefDataFrame, muRefCovariateDataFrame),
                 quote(ka <- exp(tka + 1)))

    expect_equal(.nm$.saemDropParameters(quote(ka <- exp(tka + eta.ka + 1)), muRefDataFrame, muRefCovariateDataFrame),
                 quote(ka <- exp(tka + 1)))

    expect_equal(.nm$.saemDropParameters(quote(ka <- exp(tka + 1 + eta.ka)), muRefDataFrame, muRefCovariateDataFrame),
                 quote(ka <- exp(tka + 1)))

    expect_equal(.nm$.saemDropParameters(quote(ka <- exp(eta.ka + 1 + tka)), muRefDataFrame, muRefCovariateDataFrame),
                 quote(ka <- exp(1 + tka)))


    expect_equal(.nm$.saemDropParameters(quote(ka <- exp(1 + eta.ka + tka)), muRefDataFrame, muRefCovariateDataFrame),
                 quote(ka <- exp(1 + tka)))


    expect_equal(.nm$.saemDropParameters(quote(ka <- exp(1 + tka + eta.ka)), muRefDataFrame, muRefCovariateDataFrame),
                 quote(ka <- exp(1 + tka)))


    expect_equal(.nm$.saemDropParameters(quote(cl <- exp(tcl + eta.cl + log(wt / 70) * cl.wt + sex * cl.sex + age * cl.age + 3)),
                                         muRefDataFrame, muRefCovariateDataFrame),
                 quote(cl <- exp(tcl + log(wt/70) * cl.wt + 3)))


    expect_equal(.nm$.saemDropParameters(quote(cl <- exp(tcl + eta.cl + log(wt / 70) * cl.wt + cl.sex * sex + age * cl.age + 3)),
                                         muRefDataFrame, muRefCovariateDataFrame),
                 quote(cl <- exp(tcl + log(wt/70) * cl.wt + 3)))

    muRefCov <- muRefCovariateDataFrame[!(muRefCovariateDataFrame$covariate %in% c("wt", "age")), ]

    expect_equal(.nm$.saemDropParameters(quote(cl <- exp(tcl + eta.cl + log(wt / 70) * cl.wt + cl.sex * sex + age * cl.age + 3)),
                                         muRefDataFrame, muRefCov),
                 quote(cl <- exp(tcl + log(wt/70) * cl.wt + age * cl.age + 3)))
  })

})
