nmTest({
  test_that("optim makes sense", {

    library(dplyr)

    dsn <- data.frame(i=1:1000) %>%
      mutate(time=exp(rnorm(n())),DV=rbinom(n(),1,exp(-1+time)/(1+exp(-1+time)))) %>%
      mutate(id=1)

    mod <- function() {
      ini({
        E0 <- 0.5
        Em <- 0.5
        E50 <- 2
        g <- fix(2)
      })
      model({
        v <- E0+Em*time^g/(E50^g+time^g)
        ll(bin) ~ DV * v - log(1 + exp(v))
      })
    }

    fit2 <- nlmixr(mod, dsn, est="optim")

    expect_true(inherits(fit2, "nlmixr2.optim"))

    fit3 <- fit2 %>% ini(g=unfix) %>% nlmixr2(dsn, "optim", nlmControl(covMethod="optim"))

    expect_true(inherits(fit3, "nlmixr2.optim"))

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

    fit1 <- nlmixr(one.cmt, nlmixr2data::theo_sd, est="optim", optimControl(method="L-BFGS-B"))

    expect_true(inherits(fit1, "nlmixr2.optim"))

  })

})
