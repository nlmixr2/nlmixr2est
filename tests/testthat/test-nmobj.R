nmTest({
  test_that("nmObject get tests", {

    one.cmt <- function() {
      ini({
        ## You may label each parameter with a comment
        tka <- 0.45 # Ka
        tcl <- log(c(0, 2.7, 100)) # Log Cl
        ## This works with interactive models
        ## You may also label the preceding line with label("label text")
        tv <- 3.45; label("log V")
        ## the label("Label name") works with all models
        eta.ka ~ 0.6
        eta.cl ~ 0.3
        eta.v ~ 0.1
        add.sd <- 0.7
      })
      model({
        ka <- exp(tka + eta.ka)
        cl <- exp(tcl + eta.cl)
        v <- exp(tv + eta.v)
        linCmt() ~ add(add.sd)
      })
    }

    fit <- .nlmixr(one.cmt, nlmixr2data::theo_sd, est="focei", control = foceiControlFast)

    expect_equal(fit$modelName, "one.cmt")

    one.cmt <- function() {
      ini({
        ## You may label each parameter with a comment
        tka <- 0.45 # Log Ka
        tcl <- log(c(0, 2.7, 100)) # Log Cl
        ## This works with interactive models
        ## You may also label the preceding line with label("label text")
        tv <- 3.45; label("log V")
        ## the label("Label name") works with all models
        add.sd <- 0.7
      })
      model({
        ka <- exp(tka)
        cl <- exp(tcl)
        v <- exp(tv)
        linCmt() ~ add(add.sd)
      })
    }

    fit2 <- .nlmixr(one.cmt, nlmixr2data::theo_sd, est="focei", control = foceiControlFast)

    expect_equal(fit2$modelName, "one.cmt")
  })
})

nmTest({
  test_that("nlmixr2() names a model the way rxode2() does", {
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
    mkModel <- function() one.cmt

    # a symbol keeps its name
    expect_equal(nlmixr(one.cmt)$modelName, "one.cmt")
    expect_equal(nlmixr(one.cmt)$modelName, rxode2::rxode2(one.cmt)$modelName)
    # a call is named by its text, not by the head of the call alone
    expect_equal(nlmixr(mkModel())$modelName, "mkModel()")
    expect_equal(nlmixr(mkModel())$modelName, rxode2::rxode2(mkModel())$modelName)
    # an anonymous model function is unnamed, not "function"
    expect_null(nlmixr(function() {
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
    })$modelName)
    # an rxUi keeps the name it has, and an unnamed one is named by its symbol
    expect_equal(nlmixr(rxode2::rxode2(one.cmt))$modelName, "one.cmt")
    ui <- rxode2::rxode2(function() {
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
    })
    expect_null(ui$modelName)
    expect_equal(nlmixr(ui)$modelName, "ui")
  })
})
