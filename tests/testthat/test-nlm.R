nmTest({
  test_that("nlm models convert strings to numbers", {

    mod <- function() {
      ini({
        E0 <- 0.5
        Em <- 0.5
        E50 <- 2
        g <- fix(2)
      })
      model({
        v <- E0+Em*time^g/(E50^g+time^g)+wt
        p <- expit(v)
        if (p < 0.5) {
          a <- "good"
        }  else {
          a <- "bad"
        }
        ll(bin) ~ DV * v - log(1 + exp(v)) + (a=="good")*0.5
      })
    }

    m <- mod()

    expect_error(rxode2::rxNorm(m$nlmRxModel$predOnly), NA)

  })

  test_that("nlm models add interp", {

    mod <- function() {
      ini({
        E0 <- 0.5
        Em <- 0.5
        E50 <- 2
        g <- fix(2)
      })
      model({
        v <- E0+Em*time^g/(E50^g+time^g)+wt
        p <- expit(v)
        ll(bin) ~ DV * v - log(1 + exp(v))
      })
    }

    m <- mod()

    expect_false(grepl("linear\\(wt\\)", rxode2::rxNorm(m$nlmRxModel$predOnly)))

    mod <- function() {
      ini({
        E0 <- 0.5
        Em <- 0.5
        E50 <- 2
        g <- fix(2)
      })
      model({
        linear(wt)
        v <- E0+Em*time^g/(E50^g+time^g)+wt
        p <- expit(v)
        ll(bin) ~ DV * v - log(1 + exp(v))
      })
    }

    m <- mod()

    expect_true(grepl("linear\\(wt\\)", rxode2::rxNorm(m$nlmRxModel$predOnly)))


  })

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



    fit2 <- .nlmixr(mod, dsn, est="nlm")

    expect_s3_class(fit2, "nlmixr2.nlm")

    fit2 <- .nlmixr(mod, dsn, est="bobyqa")

    expect_s3_class(fit2, "nlmixr2.bobyqa")

    fit2 <- .nlmixr(mod, dsn, est="uobyqa")

    expect_s3_class(fit2, "nlmixr2.uobyqa")

    fit2 <- .nlmixr(mod, dsn, est="newuoa")

    expect_s3_class(fit2, "nlmixr2.newuoa")

    fit2 <- .nlmixr(mod, dsn, est="n1qn1")

    expect_s3_class(fit2, "nlmixr2.n1qn1")

    fit2 <- .nlmixr(mod, dsn, est="lbfgsb3c")

    expect_s3_class(fit2, "nlmixr2.lbfgsb3c")

    fit3 <- fit2 |>
      ini(g=unfix) |>
      .nlmixr(dsn, "nlm", nlmControl(solveType="grad"))

    expect_s3_class(fit3, "nlmixr2.nlm")

    fit4 <- fit2 |>
      ini(g=unfix) |>
      .nlmixr(dsn, "nlm", nlmControl(solveType="fun"))

    expect_s3_class(fit4, "nlmixr2.nlm")

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

    fit2 <- .nlmixr(one.cmt, nlmixr2data::theo_sd, est="nlm")

    fit1 <- .nlmixr(one.cmt, nlmixr2data::theo_sd, est="nlm",
                   nlmControl(scaleTo=0.0, scaleType="multAdd"))

    expect_s3_class(fit1, "nlmixr2.nlm")
  })
})
