nmTest({
  test_that("optim makes sense", {

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
        ll(bin) ~ DV * v - log(1 + exp(v))
      })
    }

    fit2 <- .nlmixr(mod, dsn, est="optim")

    expect_true(inherits(fit2, "nlmixr2.optim"))

    fit3 <- fit2 |>
      ini(g=unfix) |>
      .nlmixr(dsn, "optim", optimControl(covMethod="optim"))

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

    fit1 <- .nlmixr(one.cmt, nlmixr2data::theo_sd, est="optim", optimControl(method="L-BFGS-B"))

    expect_true(inherits(fit1, "nlmixr2.optim"))

  })

  test_that("optim sugar aliases dispatch and honor bounds", {
    # every alias is registered as an nlmixr2Est S3 method
    for (.a in c("neldermead", "bfgs", "cg", "lbfgsb", "sann", "brent")) {
      expect_true(!is.null(utils::getS3method("nlmixr2Est", .a)),
                  info = .a)
    }
    # bounded-capable methods report unbounded=FALSE, the rest TRUE
    .unb <- function(a) attr(utils::getS3method("nlmixr2Est", a), "unbounded")(NULL)
    expect_false(.unb("brent"))
    expect_false(.unb("lbfgsb"))
    for (.a in c("neldermead", "bfgs", "cg", "sann")) {
      expect_true(.unb(.a), info = .a)
    }

    # the alias forces the corresponding optim() method end-to-end
    one.cmt <- function() {
      ini({ tka <- 0.45; tcl <- 1; tv <- 3.45; add.sd <- 0.7 })
      model({ ka <- exp(tka); cl <- exp(tcl); v <- exp(tv)
        linCmt() ~ add(add.sd) })
    }
    fitNM <- .nlmixr(one.cmt, nlmixr2data::theo_sd, est = "neldermead")
    expect_equal(fitNM$optimControl$method, "Nelder-Mead")
    fitLB <- .nlmixr(one.cmt, nlmixr2data::theo_sd, est = "lbfgsb")
    expect_equal(fitLB$optimControl$method, "L-BFGS-B")
  })

})
