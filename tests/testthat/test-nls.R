nmTest({
  test_that("nls makes sense", {

    d <- nlmixr2data::theo_sd

    d <- d[d$AMT != 0 | d$DV != 0, ]

    one.cmt <- function() {
      ini({
        tka <- 0.45
        tcl <- log(c(0, 2.7, 100))
        tv <- 3.45
        add.sd <- 0.7
      })
      model({
        ka <- exp(tka)
        cl <- exp(tcl)
        v <- exp(tv)
        linCmt() ~ add(add.sd)
      })
    }

    fit1 <- nlmixr(one.cmt, d, est="nls")

    expect_true(inherits(fit1, "nlmixr2.nls"))

    Treated <- Puromycin[Puromycin$state == "treated", ]
    names(Treated) <- gsub("rate", "DV", gsub("conc", "time", names(Treated)))
    Treated$ID <- 1

    f <- function() {
      ini({
        Vm <- 200
        K <- 0.1
        prop.sd <- 0.1
      })
      model({
        pred <- (Vm * time)/(K + time)
        pred ~ prop(prop.sd)
      })
    }

    fit1 <- nlmixr(f, Treated, est="nls", control=nlsControl(algorithm="default"))

    expect_true(inherits(fit1, "nlmixr2.nls"))
  })
})
