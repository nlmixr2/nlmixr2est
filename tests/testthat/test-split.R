.b <- loadNamespace("nlmixr2est")

test_that("test pure mu ref", {
  expect_equal(.b$.getPureMuRef(quote(cl <- tcl),
                                muRefCurEval=data.frame(parameter="tcl", curEval="",
                                                        low=NA_character_, hi=NA_character_)),
               c(tcl="cl"))

  expect_equal(.b$.getPureMuRef(quote(cl <- tcl),
                                muRefCurEval=data.frame(parameter="tcl", curEval="exp",
                                                        low=NA_character_, hi=NA_character_)), NULL)

  expect_equal(.b$.getPureMuRef(quote(cl <- exp(tcl)),
                                muRefCurEval=data.frame(parameter="tcl", curEval="exp",
                                                        low=NA_character_, hi=NA_character_)),
               c(tcl="cl"))

  expect_equal(.b$.getPureMuRef(quote(cl <- exp(tcl)),
                                muRefCurEval=data.frame(parameter="tcl", curEval="",
                                                        low=NA_character_, hi=NA_character_)),
               NULL)


  expect_equal(.b$.getPureMuRef(quote(cl <- expit(tcl, 0, 1)),
                                muRefCurEval=data.frame(parameter="tcl", curEval="expit",
                                                        low=NA_character_, hi=NA_character_)),
               c(tcl="cl"))


  expect_equal(.b$.getPureMuRef(quote(cl <- expit(tcl, 0, 2)),
                                muRefCurEval=data.frame(parameter="tcl", curEval="expit",
                                                        low=NA_character_, hi=NA_character_)),
               NULL)


  expect_equal(.b$.getPureMuRef(quote(cl <- expit(tcl, 0, 2)),
                                muRefCurEval=data.frame(parameter="tcl", curEval="expit",
                                                        low=0, hi=2)),
               c(tcl="cl"))


  expect_equal(.b$.getPureMuRef(quote(cl <- expit(tcl, 0.5, 1)),
                                muRefCurEval=data.frame(parameter="tcl", curEval="expit",
                                                        low=NA_character_, hi=NA_character_)),
               NULL)

  expect_equal(.b$.getPureMuRef(quote(cl <- expit(tcl, 0.5, 1)),
                                muRefCurEval=data.frame(parameter="tcl", curEval="expit",
                                                        low=0.5, hi=NA_character_)),
               c(tcl="cl"))


  expect_equal(.b$.getPureMuRef(quote(cl <- expit(tcl, 0.5)),
                                muRefCurEval=data.frame(parameter="tcl", curEval="expit",
                                                        low=0.5, hi=NA_character_)),
               c(tcl="cl"))


  expect_equal(.b$.getPureMuRef(quote(cl <- expit(tcl, 0.5)),
                                muRefCurEval=data.frame(parameter="tcl", curEval="expit",
                                                        low=0.4, hi=NA_character_)),
               NULL)

  expect_equal(.b$.getPureMuRef(quote(cl <- expit(tcl)),
                                muRefCurEval=data.frame(parameter="tcl", curEval="expit",
                                                        low=0, hi=1)),
               c(tcl="cl"))

  expect_equal(.b$.getPureMuRef(quote(cl <- expit(tcl)),
                                muRefCurEval=data.frame(parameter="tcl", curEval="expit",
                                                        low=NA_real_, hi=1)),
               c(tcl="cl"))

  expect_equal(.b$.getPureMuRef(quote(cl <- expit(tcl)),
                                muRefCurEval=data.frame(parameter="tcl", curEval="expit",
                                                        low=0, hi=NA_real_)),
               c(tcl="cl"))

  expect_equal(.b$.getPureMuRef(quote(cl(0) <- tcl),
                                muRefCurEval=data.frame(parameter="tcl", curEval="",
                                                        low=NA_real_, hi=NA_real_)),
               NULL)
})

test_that("test split", {
  PK_1cmt <- function() {
    description <- "One compartment PK model with linear clearance"
    ini({
      lka <- 0.45 ; label("Absorption rate (Ka)")
      lcl <- 1 ; label("Clearance (CL)")
      lvc  <- 3.45 ; label("Central volume of distribution (V)")
      prop.err <- 0.5 ; label("Proportional residual error (fraction)")
    })
    model({
      ka <- exp(lka)
      cl <- exp(lcl)
      vc  <- exp(lvc)

      cp <- linCmt()
      cp ~ prop(prop.err)
    })
  }
  
  mod <- PK_1cmt()
  expect_error(mod$getSplitMuModel, NA)
})


test_that("non mu-referenced split works correctly", {
  f <- function() {
    ini({
      tke <- 0.5
      eta.ke ~ 0.04
      prop.sd <- sqrt(0.1)
    })
    model({
      ke <- tke * exp(eta.ke)
      ipre <- 10 * exp(-ke * t)
      f2 <- ipre / (ipre + 5)
      f3 <- f2 * 3
      lipre <- log(ipre)
      ipre ~ prop(prop.sd)
    })
  }

  ui <- f()

  expect_error(ui$getSplitMuModel, NA)
})
