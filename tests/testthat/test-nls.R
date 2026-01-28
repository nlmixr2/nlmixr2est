nmTest({
  test_that("nls supports interp", {

    one.cmt <- function() {
      ini({
        tka <- fix(0.45)
        tcl <- log(c(0, 2.7, 100))
        tv <- 3.45
        add.sd <- 0.7
      })
      model({
        ka <- exp(tka)
        cl <- exp(tcl) + wt
        v <- exp(tv)
        linCmt() ~ add(add.sd)
      })
    }

    f <- one.cmt()

    expect_false(grepl("linear\\(wt\\)", rxode2::rxNorm(f$nlsRxModel$predOnly)))

    one.cmt <- function() {
      ini({
        tka <- fix(0.45)
        tcl <- log(c(0, 2.7, 100))
        tv <- 3.45
        add.sd <- 0.7
      })
      model({
        linear(wt)
        ka <- exp(tka)
        cl <- exp(tcl) + wt
        v <- exp(tv)
        linCmt() ~ add(add.sd)
      })
    }

    f <- one.cmt()

    expect_true(grepl("linear\\(wt\\)", rxode2::rxNorm(f$nlsRxModel$predOnly)))


  })

  test_that("nls all 1 issue", {

    pheno <- function() {
      ini({
        tcl <- log(1) # typical value of clearance
        tv <-  log(1)   # typical value of volume
        add.err <- 0.1    # residual variability
      })
      model({
        cl <- exp(tcl ) # individual value of clearance
        v <- exp(tv)    # individual value of volume
        ke <- cl / v            # elimination rate constant
        d/dt(A1) = - ke * A1    # model differential equation
        cp = A1 / v             # concentration in plasma
        cp ~ add(add.err)       # define error model
      })
    }

    expect_error(.nlmixr(pheno, nlmixr2data::pheno_sd, est="nls",
                         nlsControl(algorithm="LM", print=0L)), NA)

  })

  test_that("nls makes sense", {

    d <- nlmixr2data::theo_sd

    d <- d[d$AMT != 0 | d$DV != 0, ]

    one.cmt <- function() {
      ini({
        tka <- fix(0.45)
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

    fit1 <- .nlmixr(one.cmt, d, est="nls", list(print=0L))

    expect_true(inherits(fit1, "nlmixr2.nls"))

    fit1 <- .nlmixr(one.cmt, d, est="nls", nlsControl(solveType = "fun",
                                                      print=0L))

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

    fit1 <- .nlmixr(f, Treated, est="nls", control=nlsControl(algorithm="default",
                                                              print=0L))

    expect_true(inherits(fit1, "nlmixr2.nls"))
  })
})
