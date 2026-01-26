nmTest({
  skip_on_cran()
  rxode2::rxUnloadAll()
  skip_on_os("windows")

  ## start looking at transformations

  ################################################################################
  ## lognormal error structures
  ################################################################################

  .lnorm <- c(40.039, 40.039, 40.039, 40.039, 40.055, 40.055, 40.042,
              40.042, 40.039, 40.039)
  .lnormProp <- c(118.419, 118.419, 118.279, 118.311, 118.311,
                  118.311, 118.221, 118.221, 118.419, 118.419)
  .lnormPropT <- c(65.886, 65.886, 65.827, 65.832, 65.832, 65.832,
                   65.803, 65.803, 65.886, 65.886)
  .lnormPropF <- c(27.293, 27.293, 27.274, 26.643, 27.035, 27.035,
                   27.253, 27.253, 27.293, 27.293)
  .lnormPow <- c(77.87, 77.87, 77.826, 77.836, 77.836, 77.836,
                 77.803, 77.803, 77.87, 77.87)
  .lnormPowT <- c(52.616, 52.616, 52.601, 52.584, 52.587, 52.587,
                  52.587, 52.587, 52.616, 52.616)

  .lnormPowF <- c(32.834, 32.834, 32.834, 32.697, 32.784, 32.784,
                  32.817, 32.817, 32.834, 32.834)
  .lnormProp1 <- c(123.318, 123.318, 123.219, 123.24, 123.24, 123.24,
                   123.166, 123.166, 123.318, 123.318)
  .lnormPropT1 <- c(81.152, 81.152, 81.13, 81.135, 81.135, 81.135,
                    81.102, 81.102, 81.152, 81.152)
  .lnormPropF1 <- c(56.785, 56.785, 56.78, 56.774, 56.775, 56.775,
                    56.771, 56.771, 56.785, 56.785)
  .lnormPropT1 <- c(81.152, 81.152, 81.13, 81.135, 81.135, 81.135,
                    81.102, 81.102, 81.152, 81.152)
  .lnormPowF1 <- c(60.494, 60.494, 60.492, 60.489, 60.49, 60.49,
                   60.485, 60.485, 60.494, 60.494)
  .lnormPow1 <- c(89.824, 89.824, 89.803, 89.809, 89.809, 89.809,
                  89.781, 89.781, 89.824, 89.824)
  .lnormProp2 <- c(118.777, 118.777, 118.646, 118.676, 118.676,
                   118.676, 118.59, 118.59, 118.777, 118.777)
  .lnormPropT2 <- c(69.981, 69.981, 69.947, 69.951, 69.951, 69.951,
                    69.924, 69.924, 69.981, 69.981)
  .lnormPropF2 <- c(45.498, 45.498, 45.495, 45.482, 45.49, 45.49,
                    45.494, 45.494, 45.498, 45.498)
  .lnormPow2 <- c(80.244, 80.244, 80.212, 80.219, 80.219, 80.219,
                  80.191, 80.191, 80.244, 80.244)
  .lnormPowF2 <- c(48.202, 48.202, 48.2, 48.193, 48.198, 48.198,
                   48.198, 48.198, 48.202, 48.202)
  .lnormPowT2 <- c(59.837, 59.837, 59.831, 59.826, 59.827, 59.827,
                   59.819, 59.819, 59.837, 59.837)
  .lnormPowT1 <- c(72.356, 72.356, 72.351, 72.351, 72.351, 72.351,
                   72.338, 72.338, 72.356, 72.356)

  testWang2007ErrorModel("lnorm", function(f) {
    f |> model(ipre ~ lnorm(lnorm.sd)) |> ini(lnorm.sd=sqrt(0.1))
  }, .lnorm, addProp = 1)

  test_that("lognormal likelihood can be determined by data too", {
    datl <- getWang2007LogData()
    suppressWarnings({
      expect_equal(setNames(round(testWang2007ErrorModel("add lnorm", function(f) {
        f |> model(lipre ~ add(add.sd)) |> ini(add.sd=sqrt(0.1))
      },
      c(-42.106, -42.106, -42.106, -42.106, -42.09, -42.09, -42.103,
        -42.103, -42.106, -42.106),
      log=TRUE) + 2 * sum(datl$DV), 3), NULL),
      .lnorm, tolerance=5e-3)
    })
  })


  ################################################################################
  ## lnorm(NA) tests
  ################################################################################

  testWang2007ErrorModel("lnorm(NA)+prop", function(f) {
    f |> model(ipre ~ lnorm(NA) + prop(prop.sd)) |> ini(prop.sd=sqrt(0.1))
  }, .lnormProp, addProp = 1)

  testWang2007ErrorModel("lnorm(NA)+propT", function(f) {
    f |> model(ipre ~ lnorm(NA) + propT(prop.sd)) |> ini(prop.sd=sqrt(0.1))
  }, .lnormPropT, addProp = 1)

  testWang2007ErrorModel("lnorm(NA)+propF", function(f) {
    f |> model(ipre ~ lnorm(NA) + propF(prop.sd, f2)) |> ini(prop.sd=sqrt(0.1))
  }, .lnormPropF, addProp = 1)

  testWang2007ErrorModel("lnorm(NA)+pow->lnorm(NA)+prop", function(f) {
    f |> model(ipre ~ lnorm(NA) + pow(prop.sd, pw)) |> ini(prop.sd=sqrt(0.1), pw=1)
  }, .lnormProp, addProp = 1)

  testWang2007ErrorModel("lnorm(NA)+powT->lnorm(NA)+propT", function(f) {
    f |> model(ipre ~ lnorm(NA) + powT(prop.sd, pw)) |> ini(prop.sd=sqrt(0.1), pw=1)
  }, .lnormPropT, addProp = 1)

  testWang2007ErrorModel("lnorm(NA)+powF->lnorm(NA)+propF", function(f) {
    f |> model(ipre ~ lnorm(NA) + powF(prop.sd, pw, f2)) |> ini(prop.sd=sqrt(0.1), pw=1)
  }, .lnormPropF, addProp = 1)

  testWang2007ErrorModel("lnorm(NA)+pow", function(f) {
    f |> model(ipre ~ lnorm(NA) + pow(prop.sd, pw)) |> ini(prop.sd=sqrt(0.1), pw=0.5)
  }, .lnormPow, addProp = 1)

  testWang2007ErrorModel("lnorm(NA)+powT", function(f) {
    f |> model(ipre ~ lnorm(NA) + powT(prop.sd, pw)) |> ini(prop.sd=sqrt(0.1), pw=0.5)
  }, .lnormPowT, addProp = 1)

  testWang2007ErrorModel("lnorm(NA)+powF", function(f) {
    f |> model(ipre ~ lnorm(NA) + powF(prop.sd, pw, f2)) |> ini(prop.sd=sqrt(0.1), pw=0.5)
  }, .lnormPowF, addProp = 1)


  ################################################################################
  ## lnorm combined1
  ################################################################################

  testWang2007ErrorModel("lnorm+prop combined1->lnorm(NA)+prop", function(f) {
    f |> model(ipre ~ lnorm(lnorm.sd) + prop(prop.sd)) |> ini(lnorm.sd=0, prop.sd=sqrt(0.1))
  }, .lnormProp, addProp = 1)

  testWang2007ErrorModel("lnorm+propT combined1->lnorm(NA)+propT", function(f) {
    f |> model(ipre ~ lnorm(lnorm.sd) + propT(prop.sd)) |> ini(lnorm.sd=0, prop.sd=sqrt(0.1))
  }, .lnormPropT, addProp = 1)

  testWang2007ErrorModel("lnorm+propF combined1->lnorm(NA)+propF", function(f) {
    f |> model(ipre ~ lnorm(lnorm.sd) + propF(prop.sd, f2)) |> ini(lnorm.sd=0, prop.sd=sqrt(0.1))
  }, .lnormPropF, addProp = 1)

  testWang2007ErrorModel("lnorm+prop combined1->lnorm(NA)", function(f) {
    f |> model(ipre ~ lnorm(lnorm.sd) + prop(prop.sd)) |> ini(lnorm.sd=sqrt(0.1), prop.sd=0)
  }, .lnorm, addProp = 1)

  testWang2007ErrorModel("lnorm+propT combined1->lnorm(NA)", function(f) {
    f |> model(ipre ~ lnorm(lnorm.sd) + propT(prop.sd)) |> ini(lnorm.sd=sqrt(0.1), prop.sd=0)
  }, .lnorm, addProp = 1)

  testWang2007ErrorModel("lnorm+propF combined1->lnorm(NA)", function(f) {
    f |> model(ipre ~ lnorm(lnorm.sd) + propF(prop.sd, f2)) |> ini(lnorm.sd=sqrt(0.1), prop.sd=0)
  }, .lnorm, addProp = 1)

  testWang2007ErrorModel("lnorm+prop combined1", function(f) {
    f |> model(ipre ~ lnorm(lnorm.sd) + prop(prop.sd)) |> ini(lnorm.sd=sqrt(0.1), prop.sd=sqrt(0.1))
  }, .lnormProp1, addProp = 1)

  testWang2007ErrorModel("lnorm+propT combined1", function(f) {
    f |> model(ipre ~ lnorm(lnorm.sd) + propT(prop.sd)) |> ini(lnorm.sd=sqrt(0.1), prop.sd=sqrt(0.1))
  }, .lnormPropT1, addProp = 1)

  testWang2007ErrorModel("lnorm+propF combined1", function(f) {
    f |> model(ipre ~ lnorm(lnorm.sd) + propF(prop.sd, f2)) |> ini(lnorm.sd=sqrt(0.1), prop.sd=sqrt(0.1))
  }, .lnormPropF1, addProp = 1)

  testWang2007ErrorModel("lnorm+powF->lnorm+propF combined1", function(f) {
    f |> model(ipre ~ lnorm(lnorm.sd) + powF(prop.sd, pw, f2)) |> ini(lnorm.sd=sqrt(0.1), prop.sd=sqrt(0.1), pw=1)
  }, .lnormPropF1, addProp = 1)

  testWang2007ErrorModel("lnorm+pow->lnorm+prop combined1", function(f) {
    f |> model(ipre ~ lnorm(lnorm.sd) + pow(prop.sd, pw)) |> ini(lnorm.sd=sqrt(0.1), prop.sd=sqrt(0.1), pw=1)
  }, .lnormProp1, addProp = 1)

  testWang2007ErrorModel("lnorm+powT->lnorm+propT combined1", function(f) {
    f |> model(ipre ~ lnorm(lnorm.sd) + powT(prop.sd, pw)) |> ini(lnorm.sd=sqrt(0.1), prop.sd=sqrt(0.1), pw=1)
  }, .lnormPropT1, addProp = 1)

  testWang2007ErrorModel("lnorm+powF combined1", function(f) {
    f |> model(ipre ~ lnorm(lnorm.sd) + powF(prop.sd, pw, f2)) |> ini(lnorm.sd=sqrt(0.1), prop.sd=sqrt(0.1), pw=0.5)
  }, .lnormPowF1, addProp = 1)

  testWang2007ErrorModel("lnorm+pow combined1", function(f) {
    f |> model(ipre ~ lnorm(lnorm.sd) + pow(prop.sd, pw)) |> ini(lnorm.sd=sqrt(0.1), prop.sd=sqrt(0.1), pw=0.5)
  }, .lnormPow1, addProp = 1)

  testWang2007ErrorModel("lnorm+powT combined1", function(f) {
    f |> model(ipre ~ lnorm(lnorm.sd) + powT(prop.sd, pw)) |> ini(lnorm.sd=sqrt(0.1), prop.sd=sqrt(0.1), pw=0.5)
  }, .lnormPowT1, addProp = 1)

  ################################################################################
  ## lnorm combined2
  ################################################################################

  testWang2007ErrorModel("lnorm+propT combined2->lnorm(NA)+propT", function(f) {
    f |> model(ipre ~ lnorm(lnorm.sd) + propT(prop.sd)) |> ini(lnorm.sd=0, prop.sd=sqrt(0.1))
  }, .lnormPropT, addProp = 2)

  testWang2007ErrorModel("lnorm+propF combined2->lnorm(NA)+propF", function(f) {
    f |> model(ipre ~ lnorm(lnorm.sd) + propF(prop.sd, f2)) |> ini(lnorm.sd=0, prop.sd=sqrt(0.1))
  }, .lnormPropF, addProp = 2)

  testWang2007ErrorModel("lnorm+prop combined2->lnorm(NA)", function(f) {
    f |> model(ipre ~ lnorm(lnorm.sd) + prop(prop.sd)) |> ini(lnorm.sd=sqrt(0.1), prop.sd=0)
  }, .lnorm, addProp = 2)

  testWang2007ErrorModel("lnorm+propT combined2->lnorm(NA)", function(f) {
    f |> model(ipre ~ lnorm(lnorm.sd) + propT(prop.sd)) |> ini(lnorm.sd=sqrt(0.1), prop.sd=0)
  }, .lnorm, addProp = 2)

  testWang2007ErrorModel("lnorm+propF combined2->lnorm(NA)", function(f) {
    f |> model(ipre ~ lnorm(lnorm.sd) + propF(prop.sd, f2)) |> ini(lnorm.sd=sqrt(0.1), prop.sd=0)
  }, .lnorm, addProp = 2)

  testWang2007ErrorModel("lnorm+pow->lnorm+prop combined2", function(f) {
    f |> model(ipre ~ lnorm(lnorm.sd) + pow(prop.sd, pw)) |> ini(lnorm.sd=sqrt(0.1), prop.sd=sqrt(0.1), pw=1)
  }, .lnormProp2, addProp = 2)

  testWang2007ErrorModel("lnorm+prop combined2", function(f) {
    f |> model(ipre ~ lnorm(lnorm.sd) + prop(prop.sd)) |> ini(lnorm.sd=sqrt(0.1), prop.sd=sqrt(0.1))
  }, .lnormProp2 ,addProp = 2)

  testWang2007ErrorModel("lnorm+propT combined2", function(f) {
    f |> model(ipre ~ lnorm(lnorm.sd) + propT(prop.sd)) |> ini(lnorm.sd=sqrt(0.1), prop.sd=sqrt(0.1))
  }, .lnormPropT2, addProp = 2)

  testWang2007ErrorModel("lnorm+propF combined1", function(f) {
    f |> model(ipre ~ lnorm(lnorm.sd) + propF(prop.sd, f2)) |> ini(lnorm.sd=sqrt(0.1), prop.sd=sqrt(0.1))
  }, .lnormPropF2, addProp = 2)

  testWang2007ErrorModel("lnorm+powF->lnorm+propF combined2", function(f) {
    f |> model(ipre ~ lnorm(lnorm.sd) + powF(prop.sd, pw, f2)) |> ini(lnorm.sd=sqrt(0.1), prop.sd=sqrt(0.1), pw=1)
  }, .lnormPropF2, addProp = 2)

  testWang2007ErrorModel("lnorm+pow->lnorm+prop combined2", function(f) {
    f |> model(ipre ~ lnorm(lnorm.sd) + pow(prop.sd, pw)) |> ini(lnorm.sd=sqrt(0.1), prop.sd=sqrt(0.1), pw=1)
  }, .lnormProp2, addProp = 2)

  testWang2007ErrorModel("lnorm+powT->lnorm+propT combined2", function(f) {
    f |> model(ipre ~ lnorm(lnorm.sd) + powT(prop.sd, pw)) |> ini(lnorm.sd=sqrt(0.1), prop.sd=sqrt(0.1), pw=1)
  }, .lnormPropT2, addProp = 2)

  testWang2007ErrorModel("lnorm+powF combined2", function(f) {
    f |> model(ipre ~ lnorm(lnorm.sd) + powF(prop.sd, pw, f2)) |> ini(lnorm.sd=sqrt(0.1), prop.sd=sqrt(0.1), pw=0.5)
  }, .lnormPowF2, addProp = 2)

  testWang2007ErrorModel("lnorm+pow combined2", function(f) {
    f |> model(ipre ~ lnorm(lnorm.sd) + pow(prop.sd, pw)) |> ini(lnorm.sd=sqrt(0.1), prop.sd=sqrt(0.1), pw=0.5)
  }, .lnormPow2, addProp = 2)

  testWang2007ErrorModel("lnorm+powT combined2", function(f) {
    f |> model(ipre ~ lnorm(lnorm.sd) + powT(prop.sd, pw)) |> ini(lnorm.sd=sqrt(0.1), prop.sd=sqrt(0.1), pw=0.5)
  }, .lnormPowT2, addProp = 2)

  rxode2::rxUnloadAll()
})
