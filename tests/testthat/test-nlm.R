nmTest({

  test_that("nlm makes sense", {

    dsn <- data.frame(i=1:1000)
    dsn$time <- exp(rnorm(1000))
    dsn$DV <- rbinom(1000,1,exp(-1+dsn$time)/(1+exp(-1+dsn$time)))

    mod <- function() {
      ini({
        E0 <- 0.5
        Em <- 0.5
        E50 <- 2
        g <- fix(2)
      })
      model({
        v <- E0+Em*time^g/(E50^g+time^g)
        p <- expit(v)
        ll(bin) ~ DV * v - log(1 + exp(v))
      })
    }

    fit2 <- nlmixr(mod, dsn, est="nlm")

    expect_true(inherits(fit2, "nlmixr2.nlm"))

    fit2 <- nlmixr(mod, dsn, est="bobyqa")

    expect_true(inherits(fit2, "nlmixr2.bobyqa"))

    fit2 <- nlmixr(mod, dsn, est="uobyqa")

    expect_true(inherits(fit2, "nlmixr2.uobyqa"))

    fit2 <- nlmixr(mod, dsn, est="newuoa")

    expect_true(inherits(fit2, "nlmixr2.newuoa"))

    fit2 <- nlmixr(mod, dsn, est="n1qn1")

    expect_true(inherits(fit2, "nlmixr2.n1qn1"))

    fit2 <- nlmixr(mod, dsn, est="lbfgsb3c")

    expect_true(inherits(fit2, "nlmixr2.lbfgsb3c"))


   fit3 <- fit2 %>% ini(g=unfix) %>% nlmixr2(dsn, "nlm", nlmControl(solveType="grad"))

    expect_true(inherits(fit3, "nlmixr2.nlm"))

    fit4 <- fit2 %>% ini(g=unfix) %>% nlmixr2(dsn, "nlm", nlmControl(solveType="fun"))

    expect_true(inherits(fit4, "nlmixr2.nlm"))

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

    skip_if_not(rxode2parse::.linCmtSens())

    fit2 <- nlmixr(one.cmt, nlmixr2data::theo_sd, est="nlm")

    fit1 <- nlmixr(one.cmt, nlmixr2data::theo_sd, est="nlm",
                   nlmControl(scaleTo=0.0, scaleType="multAdd"))

    expect_true(inherits(fit1, "nlmixr2.nlm"))

  })

})
