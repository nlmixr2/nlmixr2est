# foceiControl(zeroTheta=) nudges a theta initialized at exactly 0 off zero
# (FOCEi scales linear thetas by |init|, which is 0 -- no scale -- at init 0).
nmTest({
  .hookB <- function(uifun, ctl = foceiControl()) {
    ui <- rxode2::rxUiDecompress(rxode2::rxode2(uifun))
    r <- .preProcessZeroTheta(ui, "focei", NULL, ctl)
    d <- (if (is.null(r)) ui else r$ui)$iniDf
    unname(d$est[d$name == "b"])
  }

  test_that("zeroTheta control round-trips and validates", {
    expect_equal(foceiControl()$zeroTheta, 0.001)
    expect_equal(foceiControl(zeroTheta = 0.01)$zeroTheta, 0.01)
    expect_error(foceiControl(zeroTheta = -1))
    expect_error(foceiControl(zeroTheta = 0))
  })

  test_that("a zero-init theta is nudged to +/- zeroTheta within bounds", {
    mUnb <- function() {
      ini({ tka<-0.45; tcl<- -1; tv<-3.45; b<-0; add.sd<-0.3; eta.cl~0.1 })
      model({ ka<-exp(tka); cl<-exp(tcl+b*WT+eta.cl); v<-exp(tv); linCmt() ~ add(add.sd) })
    }
    # unbounded -> +zeroTheta
    expect_equal(.hookB(mUnb), 0.001)

    # +zeroTheta out of bounds, -zeroTheta in -> -zeroTheta
    mNeg <- function() {
      ini({ tka<-0.45; tcl<- -1; tv<-3.45; b<-c(-1, 0, 0.0005); add.sd<-0.3; eta.cl~0.1 })
      model({ ka<-exp(tka); cl<-exp(tcl+b*WT+eta.cl); v<-exp(tv); linCmt() ~ add(add.sd) })
    }
    expect_equal(.hookB(mNeg), -0.001)

    # neither +/- zeroTheta within bounds -> error
    mBad <- function() {
      ini({ tka<-0.45; tcl<- -1; tv<-3.45; b<-c(-0.0005, 0, 0.0005); add.sd<-0.3; eta.cl~0.1 })
      model({ ka<-exp(tka); cl<-exp(tcl+b*WT+eta.cl); v<-exp(tv); linCmt() ~ add(add.sd) })
    }
    expect_error(.hookB(mBad), "off 0")

    # a theta fixed at 0 is left untouched
    mFix <- function() {
      ini({ tka<-0.45; tcl<- -1; tv<-3.45; b<-fix(0); add.sd<-0.3; eta.cl~0.1 })
      model({ ka<-exp(tka); cl<-exp(tcl+b*WT+eta.cl); v<-exp(tv); linCmt() ~ add(add.sd) })
    }
    expect_equal(.hookB(mFix), 0)
  })

  test_that("linear thetas get guarded 1/|init|, log thetas 1", {
    scB <- function(binit) {
      # bare ui has no scaleCband control -> the default c(0.1, 10) band applies
      m <- eval(bquote(function() {
        ini({ tka<-0.45; tcl<- -1; tv<-3.45; b<-.(binit); add.sd<-0.3; eta.cl~0.1 })
        model({ ka<-exp(tka); cl<-exp(tcl+b*WT+eta.cl); v<-exp(tv); linCmt() ~ add(add.sd) })
      }))
      ui <- rxode2::rxUiDecompress(rxode2::rxode2(m))
      sc <- ui$scaleCtheta
      names(sc) <- ui$iniDf$name[!is.na(ui$iniDf$ntheta) & !ui$iniDf$fix]
      unname(sc["b"])
    }
    # in-band: 1/|init| kept (results preserved).  b=2.5 -> 1/2.5 = 0.4
    expect_equal(scB(2.5), 0.4)
    # small init: 1/0.02 = 50 > 10 -> fall back to |init| = 0.02
    expect_equal(scB(0.02), 0.02)
    # large init: 1/40 = 0.025 < 0.1 -> fall back to |init| = 40
    expect_equal(scB(40), 40)
    # log-scaled structural theta -> 1
    m <- function() {
      ini({ tka<-0.45; tcl<- -1; tv<-3.45; b<-2.5; add.sd<-0.3; eta.cl~0.1 })
      model({ ka<-exp(tka); cl<-exp(tcl+b*WT+eta.cl); v<-exp(tv); linCmt() ~ add(add.sd) })
    }
    ui <- rxode2::rxUiDecompress(rxode2::rxode2(m))
    sc <- ui$scaleCtheta
    names(sc) <- ui$iniDf$name[!is.na(ui$iniDf$ntheta) & !ui$iniDf$fix]
    expect_equal(unname(sc["tka"]), 1)
  })
})
