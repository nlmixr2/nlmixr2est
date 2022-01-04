nlmixr2Test({
  test_that("theo tests", {

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

    f <- nlmixr(one.cmt)

    expect_equal(f$saemModel0,
                 quote(rxModelVars({
                   ka <- exp(tka)
                   cl <- exp(tcl)
                   v <- exp(tv)
                   rx_pred_ <- linCmt()
                 })))

    expect_equal(f$saemParamsToEstimate,
                 c("tka", "tcl", "tv"))

    expect_equal(f$saemParams,
                 "params(tka,tcl,tv)")

    expect_equal(f$saemInParsAndMuRefCovariates,
                 list(inPars = character(0), covars = character(0)))

    expect_equal(f$saemFixed,
                 c(tka = FALSE, tcl = FALSE, tv = FALSE))

    expect_equal(f$saemEtaTrans,
                 1:3)

    expect_equal(f$saemOmegaTrans,
                 1:3)

    expect_equal(f$saemModelOmega,
                 structure(c(1, 0, 0, 0, 1, 0, 0, 0, 1), .Dim = c(3L, 3L)))

    expect_equal(f$saemModelOmegaFixed,
                 matrix(rep(0, 9), 3, 3))

    expect_equal(f$saemModelOmegaFixedValues,
                 structure(c(0.6, 0, 0, 0, 0.3, 0, 0, 0, 0.1), .Dim = c(3L, 3L)))

    expect_equal(f$saemLow, -Inf)

    expect_equal(f$saemHi, -Inf)

    expect_equal(f$saemPropT, 0) # not proportional on transformed scale

    expect_equal(f$saemYj, 2) # normal translation

    expect_equal(f$saemResMod, 1) # additive

    # keep all omegas
    expect_equal(f$saemParHistOmegaKeep,
                 c(eta.ka = 1L, eta.cl = 1L, eta.v = 1L))

    expect_equal(f$saemParHistNames,
                 c("tka", "tcl", "tv", "V(eta.ka)", "V(eta.cl)", "V(eta.v)", "add.sd"))

    # y= f + (a + b*f^c)*varepsilon # combined1
    # y= f + sqrt(a^2 + b^2*f^(2*c))*varepsilon # combined2

    # Ares
    expect_equal(f$saemAres, 0.7)

    # Bres
    expect_equal(f$saemBres, 1)

    # Cres
    expect_equal(f$saemCres, 1)

    # lambda for boxCox or yeoJohnson
    expect_equal(f$saemLres, 1)

    expect_equal(f$saemLogEta, c(tka = TRUE, tcl = TRUE, tv = TRUE))

    expect_equal(f$saemInitTheta,
                 structure(c(1.56831218549017, 2.7, 31.5003923087479), .Names = c("", "", "")))

    expect_equal(f$saemInitOmega,
                 c(tka = 0.6, tcl = 0.3, tv = 0.1))

    expect_equal(f$saemThetaDataFrame,
                 structure(list(lower = c(-Inf, -Inf, -Inf, -Inf), theta = c(0.45, 0.993251773010283, 3.45, 0.7), fixed = c(FALSE, FALSE, FALSE, FALSE), upper = c(Inf, Inf, Inf, Inf)), class = "data.frame", row.names = c("tka", "tcl", "tv", "add.sd")))

  })
})
