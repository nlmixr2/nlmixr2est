nmTest({
  skip_on_cran()
  rxode2::rxUnloadAll()
  skip_on_os("windows")

  ## Test lambda=0.5

  ################################################################################
  ## BoxCox(0.5)
  ################################################################################

  .boxCoxAdd <- c(10.402, 10.402, 10.402, 10.402, 11.078, 11.078,
                  10.448, 10.448, 10.402, 10.402)

  testWang2007ErrorModel("boxCox(0.5)+add", function(f) {
    f |> model(ipre ~ add(add.sd) + boxCox(lambda)) |>
      ini(add.sd=sqrt(0.1), lambda=0.5)
  }, .boxCoxAdd, addProp = 1)

  ################################################################################
  ## boxCox(0.5)+prop
  ################################################################################

  .boxCoxProp <- c(77.857, 77.857, 77.698, 77.732, 77.732, 77.732,
                   77.679, 77.679, 77.857, 77.857)
  testWang2007ErrorModel("boxCox(0.5)+prop", function(f) {
    f |> model(ipre ~ boxCox(lambda) + prop(prop.sd)) |> ini(prop.sd=sqrt(0.1), lambda=0.5)
  }, .boxCoxProp, addProp = 1)

  .boxCoxPropT <- c(48.455, 48.455, 48.321, 48.287, 48.291, 48.291,
                    48.352, 48.352, 48.455, 48.455)
  testWang2007ErrorModel("boxCox(0.5)+propT", function(f) {
    f |> model(ipre ~ boxCox(lambda) + propT(prop.sd)) |>
      ini(prop.sd=sqrt(0.1), lambda=0.5)
  }, .boxCoxPropT, addProp = 1)

  .boxCoxPropF <- c(4.243, 4.243, 4.216, 3.537, 6.584, 6.584, 4.371,
                    4.371, 4.243, 4.243)
  testWang2007ErrorModel("boxCox(0.5)+propF", function(f) {
    f |> model(ipre ~ boxCox(lm) + propF(prop.sd, f2)) |>
      ini(prop.sd=sqrt(0.1), lm=0.5)
  }, .boxCoxPropF, addProp = 1)

  testWang2007ErrorModel("boxCox(0.5)+pow->boxCox(0.5)+prop", function(f) {
    f |> model(ipre ~ boxCox(lm) + pow(prop.sd, pw)) |>
      ini(prop.sd=sqrt(0.1), pw=1, lm=0.5)
  }, .boxCoxProp, addProp = 1)

  testWang2007ErrorModel("boxCox(0.5)+powT->boxCox(0.5)+propT", function(f) {
    f |> model(ipre ~ boxCox(lm) + powT(prop.sd, pw)) |>
      ini(prop.sd=sqrt(0.1), pw=1, lm=0.5)
  }, .boxCoxPropT, addProp = 1)

  testWang2007ErrorModel("boxCox(0.5)+powF->boxCox(0.5)+propF", function(f) {
    f |> model(ipre ~ boxCox(lm) + powF(prop.sd, pw, f2)) |>
      ini(prop.sd=sqrt(0.1), pw=1, lm=0.5)
  }, .boxCoxPropF, addProp = 1)

  .boxCoxPow <- c(39.674, 39.674, 39.629, 39.563, 39.573, 39.573,
                  39.64, 39.64, 39.674, 39.674)

  testWang2007ErrorModel("boxCox(0.5)+pow", function(f) {
    f |> model(ipre ~ boxCox(lm) + pow(prop.sd, pw)) |>
      ini(prop.sd=sqrt(0.1), pw=0.5, lm=0.5)
  }, .boxCoxPow, addProp = 1)

  .boxCoxPowT <- c(27.311, 27.311, 27.296, 27.082, 27.148, 27.148,
                   27.297, 27.297, 27.311, 27.311)
  testWang2007ErrorModel("boxCox(0.5)+powT", function(f) {
    f |> model(ipre ~ boxCox(lm) + powT(prop.sd, pw)) |>
      ini(prop.sd=sqrt(0.1), pw=0.5, lm=0.5)
  }, .boxCoxPowT, addProp = 1)

  .boxCoxPowF <- c(7.046, 7.046, 7.032, 6.617, 8.287, 8.287, 7.118,
                   7.118, 7.046, 7.046)
  testWang2007ErrorModel("boxCox(0.5)+powF", function(f) {
    f |> model(ipre ~ boxCox(lm) + powF(prop.sd, pw, f2)) |>
      ini(prop.sd=sqrt(0.1), pw=0.5, lm=0.5)
  }, .boxCoxPowF, addProp = 1)


  ################################################################################
  ## lnorm combined1
  ################################################################################

  testWang2007ErrorModel("boxCox(0.5)+add+prop combined1->boxCox(0.5)+prop", function(f) {
    f |> model(ipre ~ add(lnorm.sd) + prop(prop.sd) + boxCox(lm)) |>
      ini(lnorm.sd=0, prop.sd=sqrt(0.1), lm=0.5)
  }, .boxCoxProp, addProp = 1)

  testWang2007ErrorModel("boxCox(0.5)+add+propT combined1->boxCox(0.5)+propT", function(f) {
    f |> model(ipre ~ add(lnorm.sd) + propT(prop.sd) + boxCox(lm)) |>
      ini(lnorm.sd=0, prop.sd=sqrt(0.1), lm=0.5)
  }, .boxCoxPropT, addProp = 1)

  testWang2007ErrorModel("boxCox(0.5)+add+propF combined1->boxCox(0.5)+propF", function(f) {
    f |> model(ipre ~ add(lnorm.sd) + propF(prop.sd, f2) + boxCox(lm)) |>
      ini(lnorm.sd=0, prop.sd=sqrt(0.1), lm=0.5)
  }, .boxCoxPropF, addProp = 1)

  testWang2007ErrorModel("boxCox(0.5)+add+prop combined1->boxCox(0.5)+add", function(f) {
    f |> model(ipre ~ add(lnorm.sd) + prop(prop.sd) + boxCox(lm)) |>
      ini(lnorm.sd=sqrt(0.1), prop.sd=0, lm=0.5)
  }, .boxCoxAdd, addProp = 1)

  testWang2007ErrorModel("boxCox(0.5)+add+propT combined1->boxCox(0.5)+add", function(f) {
    f |> model(ipre ~ add(lnorm.sd) + propT(prop.sd) + boxCox(lm)) |>
      ini(lnorm.sd=sqrt(0.1), prop.sd=0, lm=0.5)
  }, .boxCoxAdd, addProp = 1)

  testWang2007ErrorModel("boxCox(0.5)+add+propF combined1->boxCox(0.5)+add", function(f) {
    f |> model(ipre ~ add(lnorm.sd) + propF(prop.sd, f2) + boxCox(lm)) |>
      ini(lnorm.sd=sqrt(0.1), prop.sd=0, lm=0.5)
  }, .boxCoxAdd, addProp = 1)

  .boxCoxAddProp1 <- c(82.622, 82.622, 82.508, 82.534, 82.534,
                       82.534, 82.481, 82.481, 82.622, 82.622)

  testWang2007ErrorModel("boxCox(0.5)+add+prop combined1", function(f) {
    f |> model(ipre ~ add(lnorm.sd) + prop(prop.sd) + boxCox(lm)) |>
      ini(lnorm.sd=sqrt(0.1), prop.sd=sqrt(0.1), lm=0.5)
  }, .boxCoxAddProp1, addProp = 1)

  .boxCoxAddPropT1 <- c(57.327, 57.327, 57.252, 57.256, 57.257,
                        57.257, 57.248, 57.248, 57.327, 57.327)

  testWang2007ErrorModel("boxCox(0.5)+add+propT combined1", function(f) {
    f |> model(ipre ~ add(lnorm.sd) + propT(prop.sd) + boxCox(lm)) |>
      ini(lnorm.sd=sqrt(0.1), prop.sd=sqrt(0.1), lm=0.5)
  }, .boxCoxAddPropT1, addProp = 1)

  .boxCoxAddPropF1 <- c(22.065, 22.065, 22.064, 21.96, 22.067,
                        22.067, 22.073, 22.073, 22.065, 22.065)

  testWang2007ErrorModel("boxCox(0.5)+add+propF combined1", function(f) {
    f |> model(ipre ~ add(lnorm.sd) + propF(prop.sd, f2) + boxCox(lm)) |>
      ini(lnorm.sd=sqrt(0.1), prop.sd=sqrt(0.1), lm=0.5)
  }, .boxCoxAddPropF1, addProp = 1)

  testWang2007ErrorModel("boxCox(0.5)+add+powF->boxCox(0.5)+add+propF combined1", function(f) {
    f |> model(ipre ~ add(lnorm.sd) + powF(prop.sd, pw, f2) + boxCox(lm)) |>
      ini(lnorm.sd=sqrt(0.1), prop.sd=sqrt(0.1), pw=1, lm=0.5)
  }, .boxCoxAddPropF1, addProp = 1)

  testWang2007ErrorModel("boxCox(0.5)+add+pow->boxCox(0.5)+add+prop combined1", function(f) {
    f |> model(ipre ~ add(lnorm.sd) + pow(prop.sd, pw) + boxCox(lm)) |>
      ini(lnorm.sd=sqrt(0.1), prop.sd=sqrt(0.1), pw=1, lm=0.5)
  }, .boxCoxAddProp1, addProp = 1)

  testWang2007ErrorModel("boxCox(0.5)+add+powT->boxCox(0.5)+add+propT combined1", function(f) {
    f |> model(ipre ~ add(lnorm.sd) + powT(prop.sd, pw) + boxCox(lm)) |>
      ini(lnorm.sd=sqrt(0.1), prop.sd=sqrt(0.1), pw=1, lm=0.5)
  }, .boxCoxAddPropT1, addProp = 1)

  .boxCoxAddPowF1 <- c(24.676, 24.676, 24.676, 24.631, 24.692,
                       24.692, 24.684, 24.684, 24.676, 24.676)

  testWang2007ErrorModel("boxCox(0.5)+add+powF combined1", function(f) {
    f |> model(ipre ~ add(lnorm.sd) + powF(prop.sd, pw, f2) + boxCox(lm)) |>
      ini(lnorm.sd=sqrt(0.1), prop.sd=sqrt(0.1), pw=0.5, lm=0.5)
  }, .boxCoxAddPowF1, addProp = 1)

  .boxCoxAddPow1 <- c(50.265, 50.265, 50.241, 50.23, 50.232, 50.232,
                      50.237, 50.237, 50.265, 50.265)

  testWang2007ErrorModel("boxCox(0.5)+add+pow combined1", function(f) {
    f |> model(ipre ~ add(lnorm.sd) + pow(prop.sd, pw) + boxCox(lm)) |>
      ini(lnorm.sd=sqrt(0.1), prop.sd=sqrt(0.1), pw=0.5, lm=0.5)
  }, .boxCoxAddPow1, addProp = 1)

  .boxCoxAddPowT1 <- c(40.45, 40.45, 40.436, 40.41, 40.416, 40.416,
                       40.434, 40.434, 40.45, 40.45)

  testWang2007ErrorModel("boxCox(0.5)+add+powT combined1", function(f) {
    f |> model(ipre ~ add(lnorm.sd) + powT(prop.sd, pw) + boxCox(lm)) |>
      ini(lnorm.sd=sqrt(0.1), prop.sd=sqrt(0.1), pw=0.5, lm=0.5)
  }, .boxCoxAddPowT1, addProp = 1)

  ################################################################################
  ## boxCox(0.5) combined2
  ################################################################################

  testWang2007ErrorModel("boxCox(0.5)+add+propT combined2->boxCox(0.5)+propT", function(f) {
    f |> model(ipre ~ add(add.sd) + propT(prop.sd) + boxCox(lambda)) |>
      ini(add.sd=0, prop.sd=sqrt(0.1), lambda=0.5)
  }, .boxCoxPropT, addProp = 2)

  testWang2007ErrorModel("boxCox(0.5)+add+propF combined2->boxCox(0.5)+propF", function(f) {
    f |> model(ipre ~ add(lnorm.sd) + boxCox(lm)+ propF(prop.sd, f2)) |>
      ini(lnorm.sd=0, prop.sd=sqrt(0.1), lm=0.5)
  }, .boxCoxPropF, addProp = 2)

  testWang2007ErrorModel("boxCox(0.5)+add+prop combined2->boxCox(0.5)+add", function(f) {
    f |> model(ipre ~ add(lnorm.sd) + prop(prop.sd) + boxCox(lm)) |>
      ini(lnorm.sd=sqrt(0.1), prop.sd=0, lm=0.5)
  }, .boxCoxAdd, addProp = 2)

  testWang2007ErrorModel("boxCox(0.5)+add+propT combined2->boxCox(0.5)+add", function(f) {
    f |> model(ipre ~ boxCox(lm) + add(lnorm.sd) + propT(prop.sd)) |>
      ini(lnorm.sd=sqrt(0.1), prop.sd=0, lm=0.5)
  }, .boxCoxAdd, addProp = 2)

  testWang2007ErrorModel("boxCox(0.5)+add+propF combined2->boxCox(0.5)+add", function(f) {
    f |> model(ipre ~ add(lnorm.sd) + boxCox(lm) + propF(prop.sd, f2)) |>
      ini(lnorm.sd=sqrt(0.1), prop.sd=0, lm=0.5)
  }, .boxCoxAdd, addProp = 2)

  .boxCoxAddProp2 <- c(78.202, 78.202, 78.052, 78.084, 78.084,
                       78.084, 78.033, 78.033, 78.202, 78.202)

  testWang2007ErrorModel("boxCox(0.5)+add+pow->boxCox(0.5)+add+prop combined2", function(f) {
    f |> model(ipre ~ boxCox(lm) + add(lnorm.sd) + pow(prop.sd, pw)) |>
      ini(lnorm.sd=sqrt(0.1), prop.sd=sqrt(0.1), pw=1, lm=0.5)
  }, .boxCoxAddProp2, addProp = 2)

  testWang2007ErrorModel("boxCox(0.5)+add+prop combined2", function(f) {
    f |> model(ipre ~ add(lnorm.sd) + prop(prop.sd) + boxCox(lm)) |>
      ini(lnorm.sd=sqrt(0.1), prop.sd=sqrt(0.1), lm=0.5)
  }, .boxCoxAddProp2 ,addProp = 2)

  .boxCoxAddPropT2 <- c(49.798, 49.798, 49.689, 49.668, 49.671,
                        49.671, 49.711, 49.711, 49.798, 49.798)

  testWang2007ErrorModel("boxCox(0.5)+add+propT combined2", function(f) {
    f |> model(ipre ~ add(lnorm.sd) + propT(prop.sd) + boxCox(lm)) |>
      ini(lnorm.sd=sqrt(0.1), prop.sd=sqrt(0.1), lm=0.5)
  }, .boxCoxAddPropT2, addProp = 2)


  .boxCoxAddPropF2 <- c(14.213, 14.213, 14.21, 14.085, 14.498,
                        14.498, 14.246, 14.246, 14.213, 14.213)

  testWang2007ErrorModel("boxCox(0.5)+add+propF combined1", function(f) {
    f |> model(ipre ~ add(lnorm.sd) + boxCox(lm)+ propF(prop.sd, f2)) |>
      ini(lnorm.sd=sqrt(0.1), prop.sd=sqrt(0.1), lm=0.5)
  }, .boxCoxAddPropF2, addProp = 2)

  testWang2007ErrorModel("boxCox(0.5)+add+powF->lnorm+propF combined2", function(f) {
    f |> model(ipre ~ add(lnorm.sd) + powF(prop.sd, pw, f2) + boxCox(lm)) |>
      ini(lnorm.sd=sqrt(0.1), prop.sd=sqrt(0.1), pw=1, lm=0.5)
  }, .boxCoxAddPropF2, addProp = 2)

  .boxCoxAddPropF2 <- c(78.202, 78.202, 78.052, 78.084, 78.084,
                        78.084, 78.033, 78.033, 78.202, 78.202)

  testWang2007ErrorModel("boxCox(0.5)+add+pow->lnorm+prop combined2", function(f) {
    f |> model(ipre ~ add(lnorm.sd) + boxCox(lm)+ pow(prop.sd, pw)) |>
      ini(lnorm.sd=sqrt(0.1), prop.sd=sqrt(0.1), pw=1, lm=0.5)
  }, .boxCoxAddPropF2, addProp = 2)

  testWang2007ErrorModel("boxCox(0.5)+add+powT->boxCox(0.5)+add+propT combined2", function(f) {
    f |> model(ipre ~ add(lnorm.sd) + powT(prop.sd, pw) + boxCox(lm)) |>
      ini(lnorm.sd=sqrt(0.1), prop.sd=sqrt(0.1), pw=1, lm=0.5)
  }, .boxCoxAddPropT2, addProp = 2)

  .boxCoxAddPowF2 <- c(15.857, 15.857, 15.856, 15.771, 16.06, 16.06,
                       15.882, 15.882, 15.857, 15.857)

  testWang2007ErrorModel("boxCox(0.5)+add+powF combined2", function(f) {
    f |> model(ipre ~ boxCox(lm) + add(lnorm.sd) + powF(prop.sd, pw, f2)) |>
      ini(lnorm.sd=sqrt(0.1), prop.sd=sqrt(0.1), pw=0.5, lm=0.5)
  }, .boxCoxAddPowF2, addProp = 2)

  .boxCoxAddPow2 <- c(41.665, 41.665, 41.631, 41.589, 41.596, 41.596,
                      41.638, 41.638, 41.665, 41.665)

  testWang2007ErrorModel("boxCox(0.5)+add+pow combined2", function(f) {
    f |> model(ipre ~ boxCox(lm) + add(lnorm.sd) + pow(prop.sd, pw)) |>
      ini(lnorm.sd=sqrt(0.1), prop.sd=sqrt(0.1), pw=0.5, lm=0.5)
  }, .boxCoxAddPow2, addProp = 2)

  .boxCoxAddPowT2 <- c(30.736, 30.736, 30.722, 30.625, 30.657,
                       30.657, 30.727, 30.727, 30.736, 30.736)

  testWang2007ErrorModel("boxCox(0.5)+add+powT combined2", function(f) {
    f |> model(ipre ~ boxCox(lm) + add(lnorm.sd) + powT(prop.sd, pw)) |>
      ini(lnorm.sd=sqrt(0.1), prop.sd=sqrt(0.1), pw=0.5, lm=0.5)
  }, .boxCoxAddPowT2, addProp = 2)

  rxode2::rxUnloadAll()
})
