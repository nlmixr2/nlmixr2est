nmTest({
  skip_on_cran()
  rxode2::rxUnloadAll()
  skip_on_os("windows")

  ## BoxCox(0) and lnorm equivalence
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
  .lnormPowF1 <- c(60.494, 60.494, 60.492, 60.489, 60.49, 60.49,
                   60.485, 60.485, 60.494, 60.494)
  .lnormPow1 <- c(89.824, 89.824, 89.803, 89.809, 89.809, 89.809,
                  89.781, 89.781, 89.824, 89.824)
  .lnormPropT1 <- c(81.152, 81.152, 81.13, 81.135, 81.135, 81.135,
                    81.102, 81.102, 81.152, 81.152)
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

  testWang2007ErrorModel("boxCox(0)+add-> lnorm", function(f) {
    f |> model(ipre ~ add(add.sd) + boxCox(lambda)) |>
      ini(add.sd=sqrt(0.1), lambda=0)
  }, .lnorm, addProp = 1)

  ################################################################################
  ## boxCox(0)+prop -> lnorm(NA) tests
  ################################################################################

  testWang2007ErrorModel("boxCox(0)+prop->lnorm(NA)+prop", function(f) {
    f |> model(ipre ~ boxCox(lambda) + prop(prop.sd)) |> ini(prop.sd=sqrt(0.1), lambda=0)
  }, .lnormProp, addProp = 1)

  testWang2007ErrorModel("boxCox(0)+propT->lnorm(NA)+propT", function(f) {
    f |> model(ipre ~ boxCox(lambda) + propT(prop.sd)) |>
      ini(prop.sd=sqrt(0.1), lambda=0)
  }, .lnormPropT, addProp = 1)

  testWang2007ErrorModel("boxCox(0)+propF->lnorm(NA)+propF", function(f) {
    f |> model(ipre ~ boxCox(lm) + propF(prop.sd, f2)) |>
      ini(prop.sd=sqrt(0.1), lm=0)
  }, .lnormPropF, addProp = 1)

  testWang2007ErrorModel("boxCox(0)+pow->lnorm(NA)+prop", function(f) {
    f |> model(ipre ~ boxCox(lm) + pow(prop.sd, pw)) |>
      ini(prop.sd=sqrt(0.1), pw=1, lm=0)
  }, .lnormProp, addProp = 1)

  testWang2007ErrorModel("boxCox(0)+powT->lnorm(NA)+propT", function(f) {
    f |> model(ipre ~ boxCox(lm) + powT(prop.sd, pw)) |>
      ini(prop.sd=sqrt(0.1), pw=1, lm=0)
  }, .lnormPropT, addProp = 1)

  testWang2007ErrorModel("boxCox(0)+powF->lnorm(NA)+propF", function(f) {
    f |> model(ipre ~ boxCox(lm) + powF(prop.sd, pw, f2)) |>
      ini(prop.sd=sqrt(0.1), pw=1, lm=0)
  }, .lnormPropF, addProp = 1)

  testWang2007ErrorModel("boxCox(0)+pow->lnorm(NA)+pow", function(f) {
    f |> model(ipre ~ boxCox(lm) + pow(prop.sd, pw)) |>
      ini(prop.sd=sqrt(0.1), pw=0.5, lm=0)
  }, .lnormPow, addProp = 1)

  testWang2007ErrorModel("boxCox(0)+powT->lnorm(NA)+powT", function(f) {
    f |> model(ipre ~ boxCox(lm) + powT(prop.sd, pw)) |>
      ini(prop.sd=sqrt(0.1), pw=0.5, lm=0)
  }, .lnormPowT, addProp = 1)

  testWang2007ErrorModel("boxCox(0)+powF->lnorm(NA)+powF", function(f) {
    f |> model(ipre ~ boxCox(lm) + powF(prop.sd, pw, f2)) |>
      ini(prop.sd=sqrt(0.1), pw=0.5, lm=0)
  }, .lnormPowF, addProp = 1)


  ################################################################################
  ## lnorm combined1
  ################################################################################

  testWang2007ErrorModel("boxCox(0)+add+prop combined1->lnorm(NA)+prop", function(f) {
    f |> model(ipre ~ add(lnorm.sd) + prop(prop.sd) + boxCox(lm)) |>
      ini(lnorm.sd=0, prop.sd=sqrt(0.1), lm=0)
  }, .lnormProp, addProp = 1)

  testWang2007ErrorModel("boxCox(0)+add+propT combined1->lnorm(NA)+propT", function(f) {
    f |> model(ipre ~ add(lnorm.sd) + propT(prop.sd) + boxCox(lm)) |>
      ini(lnorm.sd=0, prop.sd=sqrt(0.1), lm=0)
  }, .lnormPropT, addProp = 1)

  testWang2007ErrorModel("boxCox(0)+add+propF combined1->lnorm(NA)+propF", function(f) {
    f |> model(ipre ~ add(lnorm.sd) + propF(prop.sd, f2) + boxCox(lm)) |>
      ini(lnorm.sd=0, prop.sd=sqrt(0.1), lm=0)
  }, .lnormPropF, addProp = 1)

  testWang2007ErrorModel("boxCox(0)+add+prop combined1->lnorm(NA)", function(f) {
    f |> model(ipre ~ add(lnorm.sd) + prop(prop.sd) + boxCox(lm)) |>
      ini(lnorm.sd=sqrt(0.1), prop.sd=0, lm=0)
  }, .lnorm, addProp = 1)

  testWang2007ErrorModel("boxCox(0)+add+propT combined1->lnorm(NA)", function(f) {
    f |> model(ipre ~ add(lnorm.sd) + propT(prop.sd) + boxCox(lm)) |>
      ini(lnorm.sd=sqrt(0.1), prop.sd=0, lm=0)
  }, .lnorm, addProp = 1)

  testWang2007ErrorModel("boxCox(0)+add+propF combined1->lnorm(NA)", function(f) {
    f |> model(ipre ~ add(lnorm.sd) + propF(prop.sd, f2) + boxCox(lm)) |>
      ini(lnorm.sd=sqrt(0.1), prop.sd=0, lm=0)
  }, .lnorm, addProp = 1)

  testWang2007ErrorModel("boxCox(0)+add+prop combined1", function(f) {
    f |> model(ipre ~ add(lnorm.sd) + prop(prop.sd) + boxCox(lm)) |>
      ini(lnorm.sd=sqrt(0.1), prop.sd=sqrt(0.1), lm=0)
  }, .lnormProp1, addProp = 1)

  testWang2007ErrorModel("boxCox(0)+add+propT combined1", function(f) {
    f |> model(ipre ~ add(lnorm.sd) + propT(prop.sd) + boxCox(lm)) |>
      ini(lnorm.sd=sqrt(0.1), prop.sd=sqrt(0.1), lm=0)
  }, .lnormPropT1, addProp = 1)

  testWang2007ErrorModel("boxCox(0)+add+propF combined1", function(f) {
    f |> model(ipre ~ add(lnorm.sd) + propF(prop.sd, f2) + boxCox(lm)) |>
      ini(lnorm.sd=sqrt(0.1), prop.sd=sqrt(0.1), lm=0)
  }, .lnormPropF1, addProp = 1)

  testWang2007ErrorModel("boxCox(0)+add+powF->lnorm+propF combined1", function(f) {
    f |> model(ipre ~ add(lnorm.sd) + powF(prop.sd, pw, f2) + boxCox(lm)) |>
      ini(lnorm.sd=sqrt(0.1), prop.sd=sqrt(0.1), pw=1, lm=0)
  }, .lnormPropF1, addProp = 1)

  testWang2007ErrorModel("boxCox(0)+add+pow->lnorm+prop combined1", function(f) {
    f |> model(ipre ~ add(lnorm.sd) + pow(prop.sd, pw) + boxCox(lm)) |>
      ini(lnorm.sd=sqrt(0.1), prop.sd=sqrt(0.1), pw=1, lm=0)
  }, .lnormProp1, addProp = 1)

  testWang2007ErrorModel("boxCox(0)+add+powT->lnorm+propT combined1", function(f) {
    f |> model(ipre ~ add(lnorm.sd) + powT(prop.sd, pw) + boxCox(lm)) |>
      ini(lnorm.sd=sqrt(0.1), prop.sd=sqrt(0.1), pw=1, lm=0)
  }, .lnormPropT1, addProp = 1)

  testWang2007ErrorModel("boxCox(0)+add+powF combined1", function(f) {
    f |> model(ipre ~ add(lnorm.sd) + powF(prop.sd, pw, f2) + boxCox(lm)) |>
      ini(lnorm.sd=sqrt(0.1), prop.sd=sqrt(0.1), pw=0.5, lm=0)
  }, .lnormPowF1, addProp = 1)

  testWang2007ErrorModel("boxCox(0)+add+pow combined1", function(f) {
    f |> model(ipre ~ add(lnorm.sd) + pow(prop.sd, pw) + boxCox(lm)) |>
      ini(lnorm.sd=sqrt(0.1), prop.sd=sqrt(0.1), pw=0.5, lm=0)
  }, .lnormPow1, addProp = 1)

  testWang2007ErrorModel("boxCox(0)+add+powT combined1", function(f) {
    f |> model(ipre ~ add(lnorm.sd) + powT(prop.sd, pw) + boxCox(lm)) |>
      ini(lnorm.sd=sqrt(0.1), prop.sd=sqrt(0.1), pw=0.5, lm=0)
  }, .lnormPowT1, addProp = 1)

  ################################################################################
  ## lnorm combined2
  ################################################################################

  testWang2007ErrorModel("boxCox(0)+add+propT combined2->lnorm(NA)+propT", function(f) {
    f |> model(ipre ~ add(add.sd) + propT(prop.sd) + boxCox(lambda)) |>
      ini(add.sd=0, prop.sd=sqrt(0.1), lambda=0)
  }, .lnormPropT, addProp = 2)

  testWang2007ErrorModel("boxCox(0)+add+propF combined2->lnorm(NA)+propF", function(f) {
    f |> model(ipre ~ add(lnorm.sd) + boxCox(lm)+ propF(prop.sd, f2)) |>
      ini(lnorm.sd=0, prop.sd=sqrt(0.1), lm=0)
  }, .lnormPropF, addProp = 2)

  testWang2007ErrorModel("boxCox(0)+add+prop combined2->lnorm(NA)", function(f) {
    f |> model(ipre ~ add(lnorm.sd) + prop(prop.sd) + boxCox(lm)) |>
      ini(lnorm.sd=sqrt(0.1), prop.sd=0, lm=0)
  }, .lnorm, addProp = 2)

  testWang2007ErrorModel("boxCox(0)+add+propT combined2->lnorm(NA)", function(f) {
    f |> model(ipre ~ boxCox(lm) + add(lnorm.sd) + propT(prop.sd)) |>
      ini(lnorm.sd=sqrt(0.1), prop.sd=0, lm=0)
  }, .lnorm, addProp = 2)

  testWang2007ErrorModel("boxCox(0)+add+propF combined2->lnorm(NA)", function(f) {
    f |> model(ipre ~ add(lnorm.sd) + boxCox(lm) + propF(prop.sd, f2)) |>
      ini(lnorm.sd=sqrt(0.1), prop.sd=0, lm=0)
  }, .lnorm, addProp = 2)

  testWang2007ErrorModel("boxCox(0)+add+pow->lnorm+prop combined2", function(f) {
    f |> model(ipre ~ boxCox(lm) + add(lnorm.sd) + pow(prop.sd, pw)) |>
      ini(lnorm.sd=sqrt(0.1), prop.sd=sqrt(0.1), pw=1, lm=0)
  }, .lnormProp2, addProp = 2)

  testWang2007ErrorModel("boxCox(0)+add+prop combined2", function(f) {
    f |> model(ipre ~ add(lnorm.sd) + prop(prop.sd) + boxCox(lm)) |>
      ini(lnorm.sd=sqrt(0.1), prop.sd=sqrt(0.1), lm=0)
  }, .lnormProp2 ,addProp = 2)

  testWang2007ErrorModel("boxCox(0)+add+propT combined2", function(f) {
    f |> model(ipre ~ add(lnorm.sd) + propT(prop.sd) + boxCox(lm)) |>
      ini(lnorm.sd=sqrt(0.1), prop.sd=sqrt(0.1), lm=0)
  }, .lnormPropT2, addProp = 2)

  testWang2007ErrorModel("boxCox(0)+add+propF combined1", function(f) {
    f |> model(ipre ~ add(lnorm.sd) + boxCox(lm)+ propF(prop.sd, f2)) |>
      ini(lnorm.sd=sqrt(0.1), prop.sd=sqrt(0.1), lm=0)
  }, .lnormPropF2, addProp = 2)

  testWang2007ErrorModel("boxCox(0)+add+powF->lnorm+propF combined2", function(f) {
    f |> model(ipre ~ add(lnorm.sd) + powF(prop.sd, pw, f2) + boxCox(lm)) |>
      ini(lnorm.sd=sqrt(0.1), prop.sd=sqrt(0.1), pw=1, lm=0)
  }, .lnormPropF2, addProp = 2)


  testWang2007ErrorModel("boxCox(0)+add+pow->lnorm+prop combined2", function(f) {
    f |> model(ipre ~ add(lnorm.sd) + boxCox(lm)+ pow(prop.sd, pw)) |>
      ini(lnorm.sd=sqrt(0.1), prop.sd=sqrt(0.1), pw=1, lm=0)
  }, .lnormProp2, addProp = 2)

  testWang2007ErrorModel("boxCox(0)+add+powT->lnorm+propT combined2", function(f) {
    f |> model(ipre ~ add(lnorm.sd) + powT(prop.sd, pw) + boxCox(lm)) |>
      ini(lnorm.sd=sqrt(0.1), prop.sd=sqrt(0.1), pw=1, lm=0)
  }, .lnormPropT2, addProp = 2)

  testWang2007ErrorModel("boxCox(0)+add+powF combined2", function(f) {
    f |> model(ipre ~ boxCox(lm) + add(lnorm.sd) + powF(prop.sd, pw, f2)) |>
      ini(lnorm.sd=sqrt(0.1), prop.sd=sqrt(0.1), pw=0.5, lm=0)
  }, .lnormPowF2, addProp = 2)

  testWang2007ErrorModel("boxCox(0)+add+pow combined2", function(f) {
    f |> model(ipre ~ boxCox(lm) + add(lnorm.sd) + pow(prop.sd, pw)) |>
      ini(lnorm.sd=sqrt(0.1), prop.sd=sqrt(0.1), pw=0.5, lm=0)
  }, .lnormPow2, addProp = 2)

  testWang2007ErrorModel("boxCox(0)+add+powT combined2", function(f) {
    f |> model(ipre ~ boxCox(lm) + add(lnorm.sd) + powT(prop.sd, pw)) |>
      ini(lnorm.sd=sqrt(0.1), prop.sd=sqrt(0.1), pw=0.5, lm=0)
  }, .lnormPowT2, addProp = 2)

  rxode2::rxUnloadAll()
})
