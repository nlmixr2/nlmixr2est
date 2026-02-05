nmTest({
  skip_on_cran()
  rxode2::rxUnloadAll()
  skip_on_os("windows")

  #################################################################################
  # Box-Cox Regression
  #################################################################################

  .propVals <- c(39.458, 39.458, 39.275, 39.207, 39.213, 39.213,
                 39.365, 39.365, 39.458, 39.458)
  .propModVals <- c(63.353, 63.353, 63.001, 63.063, 63.063, 63.063,
                    63.063, 63.063, 63.353, 63.353)
  .propFVals <- c(6.496, 6.496, 6.488, 6.275, 9.262, 9.262, 6.545,
                  6.545, 6.496, 6.496)
  .propFModVals <- c(19.177, 19.177, 19.07, 18.202, 18.333, 18.333,
                     19.158, 19.158, 19.177, 19.177)
  .addVals <- c(-2.059, -2.059, -2.059, -2.059, 0.026, 0.026, -1.997,
                -1.997, -2.059, -2.059)
  .addModVals <- c(3.238, 3.238, 3.207, 2.438, 3.311, 3.311, 3.298,
                   3.298, 3.238, 3.238)
  .powVals <- c(9.966, 9.966, 9.948, 9.331, 9.651, 9.651, 10.007,
                10.007, 9.966, 9.966)
  .powF1Vals <- c(27.301, 27.301, 27.147, 26.854, 26.888, 26.888,
                  27.259, 27.259, 27.301, 27.301)
  .powF2Vals <- c(17.831, 17.831, 17.785, 16.617, 16.852, 16.852,
                  17.871, 17.871, 17.831, 17.831)
  .powF3Vals <- c(79.733, 79.733, 79.371, 79.448, 79.448, 79.448,
                  79.389, 79.389, 79.733, 79.733)
  .powFModVals <- c(0.776, 0.776, 0.772, 0.58, 3.2, 3.2, 0.828,
                    0.828, 0.776, 0.776)
  .powFMod1Vals <- c(24.877, 24.877, 24.773, 24.518, 24.554, 24.554,
                     24.842, 24.842, 24.877, 24.877)
  .powFMod2Vals <- c(49.312, 49.312, 49.311, 49.296, 52.603, 52.603,
                     49.34, 49.34, 49.312, 49.312)
  .powFMod3Vals <- c(10.848, 10.848, 10.784, 9.446, 9.96, 9.96,
                     10.915, 10.915, 10.848, 10.848)
  .boxCoxPropTVals <- c(35.132, 35.132, 34.879, 34.694, 34.709,
                        34.709, 35.036, 35.036, 35.132, 35.132)
  .boxCoxPropTModVals <- c(58.21, 58.21, 57.741, 57.802, 57.803,
                           57.803, 57.884, 57.884, 58.21, 58.21)
  .boxCoxPropTVals <- c(35.132, 35.132, 34.879, 34.694, 34.709,
                        34.709, 35.036, 35.036, 35.132, 35.132)
  .boxCoxPowTVals <- c(9.11, 9.11, 9.095, 8.215, 8.644, 8.644,
                       9.152, 9.152, 9.11, 9.11)
  .boxCoxPowTF1Vals <- c(25.57, 25.57, 25.388, 24.951, 25.002,
                         25.002, 25.53, 25.53, 25.57, 25.57)
  .boxCoxPowTF2Vals <- c(16.507, 16.507, 16.423, 14.981, 15.314,
                         15.314, 16.568, 16.568, 16.507, 16.507)
  .boxCoxPowTF3Vals <- c(76.496, 76.496, 76.097, 76.179, 76.179,
                         76.179, 76.137, 76.137, 76.496, 76.496)
  .boxCoxAddModPropTMod2 <- c(100.737, 100.737, 100.351, 100.437,
                              100.437, 100.437, 100.293, 100.293, 100.737, 100.737)
  .boxCoxAddPropT2Vals <- c(35.457, 35.457, 35.224, 35.056, 35.07,
                            35.07, 35.367, 35.367, 35.195, 35.195)
  .boxCoxPropT2ModVals <- c(58.21, 58.21, 57.741, 57.802, 57.803,
                            57.803, 57.884, 57.884, 58.21, 58.21)
  .boxCoxPropAddPropT2Vals <- c(100.737, 100.737, 100.351, 100.437, 100.437, 100.437)
  .addProp2 <- c(39.735, 39.735, 39.562, 39.499, 39.505, 39.505,
                 39.647, 39.647, 39.735, 39.735)
  .addModPropMod2 <- c(106.308, 106.308, 106.013, 106.079, 106.079,
                       106.079, 105.948, 105.948, 106.308, 106.308)
  .addModPropFModVals2 <- c(54.317, 54.317, 54.14, 54.165, 54.166,
                            54.166, 54.148, 54.148, 54.317, 54.317)
  .addPropFVals2 <- c(-2.321, -2.321, -2.322, -2.454, -0.65, -0.65,
                      -2.247, -2.247, -2.321, -2.321)
  .addModPropMod1 <- c(106.308, 106.308, 106.013, 106.079, 106.079,
                       106.079, 105.948, 105.948, 106.308, 106.308)
  .addProp1 <- c(43.554, 43.554, 43.416, 43.394, 43.398, 43.398,
                 43.469, 43.469, 43.554, 43.554)

  testWang2007ErrorModel("prop+boxCox->prop", function(f) {
    f |> model(ipre ~ prop(prop.sd) + boxCox(lambda)) |> ini(prop.sd=sqrt(0.1), lambda=1)
  }, .propVals)

  testWang2007ErrorModel("prop+boxCox->propMod", function(f) {
    f |> model(ipre ~ prop(f2) + boxCox(lambda)) |> ini(lambda=1)
  }, .propModVals)

  # In this case propT = prop
  testWang2007ErrorModel("propT+boxCox->propT", function(f) {
    f |> model(ipre ~ propT(prop.sd) + boxCox(lambda)) |> ini(prop.sd=sqrt(0.1), lambda=1)
  }, .boxCoxPropTVals)

  testWang2007ErrorModel("propTMod+boxCox->propTMod", function(f) {
    f |> model(ipre ~ propT(f2) + boxCox(lambda)) |> ini(lambda=1)
  }, .boxCoxPropTModVals)

  testWang2007ErrorModel("propF+boxCox->propF", function(f) {
    f |> model(ipre ~ propF(prop.sd, f2) + boxCox(lambda)) |> ini(prop.sd=sqrt(0.1), lambda=1)
  }, .propFVals)

  testWang2007ErrorModel("propFMod+boxCox->propFMod", function(f) {
    f |> model(ipre ~ propF(lipre, f2) + boxCox(lambda)) |> ini(lambda=1)
  }, .propFModVals)

  ################################################################################
  # boxCox -> Additive Model Tests
  ################################################################################
  testWang2007ErrorModel("add+boxCox->add", function(f) {
    f |> model(ipre ~ add(add.sd) + boxCox(lambda)) |> ini(add.sd=sqrt(0.1), lambda=1)
  }, .addVals)

  testWang2007ErrorModel("addMod+boxCox->addMod", function(f) {
    f |> model(ipre ~ add(f2) + boxCox(lambda)) |> ini(lambda=1)
  }, .addModVals)

  ################################################################################
  # boxCox -> Power Model Tests
  ################################################################################
  testWang2007ErrorModel("boxCox+pow1->prop", function(f) {
    f |> model(ipre ~ pow(pow.sd, pw) + boxCox(lambda)) |> ini(pow.sd=sqrt(0.1), pw=1, lambda=1)
  }, .propVals)

  testWang2007ErrorModel("boxCox+powT1=propT=prop", function(f) {
    f |> model(ipre ~ powT(pow.sd, pw) + boxCox(lambda)) |> ini(pow.sd=sqrt(0.1), pw=1, lambda=1)
  }, .boxCoxPropTVals)

  testWang2007ErrorModel("boxCox+powF1=propF", function(f) {
    f |> model(ipre ~ powF(pow.sd, pw, f2) + boxCox(lambda)) |> ini(pow.sd=sqrt(0.1), pw=1, lambda=1)
  }, .propFVals)

  testWang2007ErrorModel("boxCox+pow1Mod=propMod", function(f) {
    f |> model(ipre ~ pow(f2, pw) + boxCox(lambda)) |> ini(pw=1, lambda=1)
  }, .propModVals)

  testWang2007ErrorModel("boxCox+pow->pow", function(f) {
    f |> model(ipre ~ pow(pow.sd, pw) + boxCox(lambda)) |> ini(pow.sd=sqrt(0.1), pw=0.5, lambda=1)
  }, .powVals)

  testWang2007ErrorModel("boxCox+powF1->powF1", function(f) {
    f |> model(ipre ~ pow(f2, pw) + boxCox(lambda)) |> ini(pw=0.5, lambda=1)
  }, .powF1Vals)

  testWang2007ErrorModel("boxCox+powF2->powF2", function(f) {
    f |> model(ipre ~ pow(pow.sd, f2) + boxCox(lambda)) |> ini(pow.sd=sqrt(0.1), lambda=1)
  }, .powF2Vals)

  testWang2007ErrorModel("boxCox+powF3->powF3", function(f) {
    f |> model(ipre ~ pow(lipre, f2) + boxCox(lambda)) |> ini(lambda=1)
  }, .powF3Vals)

  testWang2007ErrorModel("boxCox+powT->powT", function(f) {
    f |> model(ipre ~ powT(pow.sd, pw) + boxCox(lambda)) |> ini(pow.sd=sqrt(0.1), pw=0.5, lambda=1)
  }, .boxCoxPowTVals)

  testWang2007ErrorModel("boxCox+powTF1->powTF1", function(f) {
    f |> model(ipre ~ powT(f2, pw) + boxCox(lambda)) |> ini(pw=0.5, lambda=1)
  }, .boxCoxPowTF1Vals)

  testWang2007ErrorModel("boxCox+powFT2->powFT2", function(f) {
    f |> model(ipre ~ powT(pow.sd, f2) + boxCox(lambda)) |> ini(pow.sd=sqrt(0.1), lambda=1)
  }, .boxCoxPowTF2Vals)

  testWang2007ErrorModel("boxCox+powTF3->powTF3", function(f) {
    f |> model(ipre ~ powT(lipre, f2) + boxCox(lambda)) |> ini(lambda=1)
  }, .boxCoxPowTF3Vals)

  testWang2007ErrorModel("boxCox+powFMod->powFMod", function(f) {
    f |> model(ipre ~ powF(pow.sd, pw, f2) + boxCox(lambda)) |> ini(pow.sd=sqrt(0.1), pw=0.5, lambda=1)
  }, .powFModVals)

  testWang2007ErrorModel("boxCox+powFMod1->powFMod1", function(f) {
    f |> model(ipre ~ powF(lipre, pw, f2) + boxCox(lambda)) |> ini(pw=0.5, lambda=1)
  }, .powFMod1Vals)

  testWang2007ErrorModel("boxCox+powFMod2->powFMod2", function(f) {
    f |> model(ipre ~ powF(pow.sd, lipre, f2) + boxCox(lambda)) |> ini(pow.sd=sqrt(0.1), lambda=1)
  }, .powFMod2Vals)

  testWang2007ErrorModel("boxCox+powFMod3->powFMod3", function(f) {
    f |> model(ipre ~ powF(lipre, f3, f2) + boxCox(lambda)) |> ini(lambda=1)
  }, .powFMod3Vals)

  ################################################################################
  # Box-Cox Add+Proportional tests (combined 2)
  ################################################################################
  testWang2007ErrorModel("boxCox+add+prop, combined 2->add", function(f) {
    f |> model(ipre ~ add(add.sd) + prop(prop.sd) + boxCox(lambda)) |> ini(add.sd=sqrt(0.1), prop.sd=0, lambda=1)
  }, .addVals, addProp = 2)

  testWang2007ErrorModel("boxCox+addMod+prop, combined 2->addMod", function(f) {
    f |> model(ipre ~ add(f2) + prop(prop.sd) + boxCox(lambda)) |> ini(prop.sd=0, lambda=1)
  }, .addModVals, addProp = 2)

  testWang2007ErrorModel("boxCox+add+prop, combined 2->prop", function(f) {
    f |> model(ipre ~ add(add.sd) + prop(prop.sd) + boxCox(lambda)) |> ini(add.sd=0, prop.sd=sqrt(0.1), lambda=1)
  }, .propVals, addProp = 2)

  testWang2007ErrorModel("boxCox+add+propMod, combined 2->propMod", function(f) {
    f |> model(ipre ~ add(add.sd) + prop(f2) + boxCox(lambda)) |> ini(add.sd=0, lambda=1)
  }, .propModVals, addProp = 2)

  testWang2007ErrorModel("boxCox+addMod+propMod, combined 2->propMod", function(f) {
    f |> model(ipre ~ add(f2) + prop(f3) + boxCox(lambda)) |> ini(lambda=1)
  }, .addModPropMod2, addProp = 2)

  testWang2007ErrorModel("boxCox+add+prop, combined 2->add+prop, combined 2", function(f) {
    f |> model(ipre ~ add(add.sd) + prop(prop.sd) + boxCox(lambda)) |> ini(add.sd=sqrt(0.1), prop.sd=sqrt(0.1), lambda=1)
  }, .addProp2, addProp = 2)

  testWang2007ErrorModel("boxCox+add+prop, combined 2->add+prop, combined 2", function(f) {
    f |> model(ipre ~ add(add.sd) + prop(prop.sd) + boxCox(lambda) + combined2()) |>
      ini(add.sd=sqrt(0.1), prop.sd=sqrt(0.1), lambda=1)
  }, .addProp2, addProp = 1)

  ## propT
  testWang2007ErrorModel("boxCox+add+propT, combined 2->add", function(f) {
    f |> model(ipre ~ add(add.sd) + propT(prop.sd) + boxCox(lambda)) |>
      ini(add.sd=sqrt(0.1), prop.sd=0, lambda=1)
  }, .addVals, addProp = 2)

  testWang2007ErrorModel("boxCox+addMod+propT, combined 2->addMod", function(f) {
    f |> model(ipre ~ add(f2) + propT(prop.sd) + boxCox(lambda)) |>
      ini(prop.sd=0, lambda=1)
  }, .addModVals, addProp = 2)

  testWang2007ErrorModel("boxCox+add+propT, combined 2->propT", function(f) {
    f |> model(ipre ~ add(add.sd) + propT(prop.sd) + boxCox(lambda)) |>
      ini(add.sd=0, prop.sd=sqrt(0.1), lambda=1)
  }, c(35.132, 35.132, 34.879, 34.694, 34.709, 34.709, 35.036, 35.036,
       35.132, 35.132),
  addProp = 2)

  testWang2007ErrorModel("boxCox+add+propTMod, combined 2->propMod", function(f) {
    f |> model(ipre ~ add(add.sd) + propT(f2) + boxCox(lambda)) |>
      ini(add.sd=0, lambda=1)
  }, .boxCoxPropTModVals, addProp = 2)

  testWang2007ErrorModel("boxCox+addMod+propTMod, combined 2->propTMod", function(f) {
    f |> model(ipre ~ add(f2) + propT(f3) + boxCox(lambda)) |>
      ini(lambda=1)
  }, .boxCoxAddModPropTMod2, addProp = 2)

  .boxCoxAddPropTc2Vals <- c(35.457, 35.457, 35.224, 35.056, 35.07,
                             35.07, 35.367, 35.367, 35.457, 35.457)

  testWang2007ErrorModel("boxCox+add+propT, combined 2", function(f) {
    f |> model(ipre ~ add(add.sd) + propT(prop.sd) + boxCox(lambda)) |>
      ini(add.sd=sqrt(0.1), prop.sd=sqrt(0.1), lambda=1)
  }, .boxCoxAddPropTc2Vals, addProp = 2)

  testWang2007ErrorModel("boxCox+add+propT, combined 2 (specified)", function(f) {
    f |> model(ipre ~ add(add.sd) + propT(prop.sd) + boxCox(lambda) + combined2()) |>
      ini(add.sd=sqrt(0.1), prop.sd=sqrt(0.1), lambda=1)
  }, .boxCoxAddPropTc2Vals, addProp = 1)

  # propF
  testWang2007ErrorModel("boxCox+add+propF, combined 2->add", function(f) {
    f |> model(ipre ~ add(add.sd) + propF(prop.sd, f2) + boxCox(lambda)) |> ini(add.sd=sqrt(0.1), prop.sd=0, lambda=1)
  }, .addVals, addProp = 2)

  testWang2007ErrorModel("boxCox+addMod+propF, combined 2->addMod", function(f) {
    f |> model(ipre ~ add(f2) + propF(prop.sd, f2) + boxCox(lambda)) |> ini(prop.sd=0, lambda=1)
  }, .addModVals, addProp = 2)

  testWang2007ErrorModel("boxCox+add+propFMod, combined 2->propMod", function(f) {
    f |> model(ipre ~ add(add.sd) + propF(lipre, f2) + boxCox(lambda)) |> ini(add.sd=0, lambda=1)
  }, .propFModVals, addProp = 2)

  testWang2007ErrorModel("boxCox+addMod+propFMod, combined 2->propMod", function(f) {
    f |> model(ipre ~ add(f2) + propF(lipre, f3) + boxCox(lambda)) |>
      ini(lambda=1)
  }, .addModPropFModVals2, addProp = 2)

  testWang2007ErrorModel("boxCox+add+propF, combined 2", function(f) {
    f |> model(ipre ~ add(add.sd) + propF(prop.sd, f2) + boxCox(lambda)) |>
      ini(add.sd=sqrt(0.1), prop.sd=sqrt(0.1), lambda=1)
  },  .addPropFVals2, addProp = 2)

  testWang2007ErrorModel("boxCox+add+propF, combined 2 (specified)", function(f) {
    f |> model(ipre ~ add(add.sd) + propF(prop.sd, f2) + boxCox(lambda) + combined2()) |>
      ini(add.sd=sqrt(0.1), prop.sd=sqrt(0.1), lambda=1)
  }, .addPropFVals2, addProp = 1)

  ################################################################################
  # Box-Box Add+Proportional tests (combined 1)
  ################################################################################
  testWang2007ErrorModel("boxCox+add+prop, combined 1->add", function(f) {
    f |> model(ipre ~ add(add.sd) + prop(prop.sd) + boxCox(lambda)) |>
      ini(add.sd=sqrt(0.1), prop.sd=0, lambda=1)
  }, .addVals, addProp = 1)

  testWang2007ErrorModel("boxCox+addMod+prop, combined 1->addMod", function(f) {
    f |> model(ipre ~ add(f2) + prop(prop.sd) + boxCox(lambda)) |>
      ini(prop.sd=0, lambda=1)
  }, .addModVals, addProp = 1)

  testWang2007ErrorModel("boxCox+add+prop, combined 1->prop", function(f) {
    f |> model(ipre ~ add(add.sd) + prop(prop.sd) + boxCox(lambda)) |>
      ini(add.sd=0, prop.sd=sqrt(0.1), lambda=1)
  }, .propVals, addProp = 1)

  testWang2007ErrorModel("boxCox+add+propMod, combined 1->propMod", function(f) {
    f |> model(ipre ~ add(add.sd) + prop(f2) + boxCox(lambda)) |>
      ini(add.sd=0, lambda=1)
  }, .propModVals, addProp = 1)

  testWang2007ErrorModel("boxCox+addMod+propMod, combined 1->propMod", function(f) {
    f |> model(ipre ~ add(f2) + prop(f3) + boxCox(lambda)) |>
      ini(lambda=1)
  },  .addModPropMod1, addProp = 2)

  testWang2007ErrorModel("boxCox+add+prop, combined 2->add+prop, combined 2", function(f) {
    f |> model(ipre ~ add(add.sd) + prop(prop.sd) + boxCox(lambda)) |>
      ini(add.sd=sqrt(0.1), prop.sd=sqrt(0.1), lambda=1)
  }, .addProp2, addProp = 2)

  testWang2007ErrorModel("boxCox+add+prop, combined 2", function(f) {
    f |> model(ipre ~ add(add.sd) + prop(prop.sd) + boxCox(lambda) + combined2()) |>
      ini(add.sd=sqrt(0.1), prop.sd=sqrt(0.1), lambda=1)
  }, .addProp2, addProp = 1)

  ## propT
  testWang2007ErrorModel("boxCox+add+propT, combined 2->add", function(f) {
    f |> model(ipre ~ add(add.sd) + propT(prop.sd) + boxCox(lambda)) |>
      ini(add.sd=sqrt(0.1), prop.sd=0, lambda=1)
  }, .addVals, addProp = 2)

  testWang2007ErrorModel("boxCox+addMod+propT, combined 2->addMod", function(f) {
    f |> model(ipre ~ add(f2) + propT(prop.sd) + boxCox(lambda)) |>
      ini(prop.sd=0, lambda=1)
  }, .addModVals, addProp = 2)

  .boxCoxPropT <- c(35.132, 35.132, 34.879, 34.694, 34.709, 34.709,
                    35.036, 35.036, 35.132, 35.132)

  testWang2007ErrorModel("boxCox+add+propT, combined 2->prop", function(f) {
    f |> model(ipre ~ add(add.sd) + propT(prop.sd) + boxCox(lambda)) |>
      ini(add.sd=0, prop.sd=sqrt(0.1), lambda=1)
  }, .boxCoxPropT, addProp = 2)

  testWang2007ErrorModel("boxCox+add+propTMod, combined 2->propMod", function(f) {
    f |> model(ipre ~ add(add.sd) + propT(f2) + boxCox(lambda)) |>
      ini(add.sd=0, lambda=1)
  }, .boxCoxPropT2ModVals, addProp = 2)
  .boxCoxAddPropT3Vals <- c(100.737, 100.737, 100.351, 100.437,
                            100.437, 100.437, 100.293, 100.293, 100.737, 100.737)

  testWang2007ErrorModel("boxCox+addMod+propTMod, combined 2->propMod", function(f) {
    f |> model(ipre ~ add(f2) + propT(f3) + boxCox(lambda)) |>
      ini(lambda=1)
  },  .boxCoxAddPropT3Vals, addProp = 2)

  .boxCoxAddPropT2Vals <- c(35.457, 35.457, 35.224, 35.056, 35.07,
                            35.07, 35.367, 35.367, 35.457, 35.457)

  testWang2007ErrorModel("boxCox+add+propT, combined 2 (specified)", function(f) {
    f |> model(ipre ~ add(add.sd) + propT(prop.sd) + boxCox(lambda) + combined2()) |>
      ini(add.sd=sqrt(0.1), prop.sd=sqrt(0.1), lambda=1)
  }, .boxCoxAddPropT2Vals, addProp = 1)

  # propF
  testWang2007ErrorModel("boxCox+add+propF, combined 2->add", function(f) {
    f |> model(ipre ~ add(add.sd) + propF(prop.sd, f2) + boxCox(lambda)) |>
      ini(add.sd=sqrt(0.1), prop.sd=0, lambda=1)
  }, .addVals, addProp = 2)

  testWang2007ErrorModel("boxCox+addMod+propF, combined 2->addMod", function(f) {
    f |> model(ipre ~ add(f2) + propF(prop.sd, f2) + boxCox(lambda)) |>
      ini(prop.sd=0, lambda=1)
  }, .addModVals, addProp = 2)

  testWang2007ErrorModel("boxCox+add+propF, combined 2->prop", function(f) {
    f |> model(ipre ~ add(add.sd) + propF(prop.sd, f2) + boxCox(lambda)) |>
      ini(add.sd=0, prop.sd=sqrt(0.1), lambda=1)
  }, .propFVals, addProp = 2)

  testWang2007ErrorModel("boxCox+add+propFMod, combined 2->propMod", function(f) {
    f |> model(ipre ~ add(add.sd) + propF(lipre, f2) + boxCox(lambda)) |>
      ini(add.sd=0, lambda=1)
  }, .propFModVals, addProp = 2)

  testWang2007ErrorModel("boxCox+addMod+propFMod, combined 2->propMod", function(f) {
    f |> model(ipre ~ add(f2) + propF(lipre, f3) + boxCox(lambda)) |>
      ini(lambda=1)
  }, .addModPropFModVals2, addProp = 2)

  testWang2007ErrorModel("boxCox+add+propF, combined 2", function(f) {
    f |> model(ipre ~ add(add.sd) + propF(prop.sd, f2) + boxCox(lambda)) |>
      ini(add.sd=sqrt(0.1), prop.sd=sqrt(0.1), lambda=1)
  },  .addPropFVals2, addProp = 2)

  testWang2007ErrorModel("boxCox+add+propF, combined 2 (specified)", function(f) {
    f |> model(ipre ~ add(add.sd) + propF(prop.sd, f2) + boxCox(lambda) + combined2()) |>
      ini(add.sd=sqrt(0.1), prop.sd=sqrt(0.1), lambda=1)
  }, .addPropFVals2, addProp = 1)

  #################################################################################################
  testWang2007ErrorModel("boxCox+add+prop, combined 1", function(f) {
    f |> model(ipre ~ add(add.sd) + prop(prop.sd) + boxCox(lambda)) |>
      ini(add.sd=sqrt(0.1), prop.sd=sqrt(0.1), lambda=1)
  }, .addProp1, addProp = 1)

  testWang2007ErrorModel("boxCox+add+prop, combined 1->add", function(f) {
    f |> model(ipre ~ add(add.sd) + prop(prop.sd) + boxCox(lambda)) |>
      ini(add.sd=sqrt(0.1), prop.sd=0, lambda=1)
  }, .addVals, addProp = 1)

  testWang2007ErrorModel("boxCox+add+prop, combined 1->prop", function(f) {
    f |> model(ipre ~ add(add.sd) + prop(prop.sd) + boxCox(lambda)) |>
      ini(add.sd=0, prop.sd=sqrt(0.1), lambda=1)
  }, .propVals, addProp = 1)

  testWang2007ErrorModel("boxCox+add+prop, combined 1 (override)", function(f) {
    f |> model(ipre ~ add(add.sd) + prop(prop.sd) + boxCox(lambda) + combined1()) |>
      ini(add.sd=sqrt(0.1), prop.sd=sqrt(0.1), lambda=1)
  }, .addProp1, addProp = 2)

  testWang2007ErrorModel("boxCox+add+prop, combined 2 (override)", function(f) {
    f |> model(ipre ~ add(add.sd) + prop(prop.sd) + boxCox(lambda) + combined2()) |>
      ini(add.sd=sqrt(0.1), prop.sd=sqrt(0.1), lambda=1)
  }, .addProp2, addProp = 1)

  ################################################################################
  # Add+Pow tests (combined 2)
  ################################################################################
  testWang2007ErrorModel("boxCox+add+pow combined 2 -> add+prop combined2", function(f) {
    f |> model(ipre ~ add(add.sd) + pow(prop.sd, pw) + boxCox(lambda)) |>
      ini(add.sd=sqrt(0.1), prop.sd=sqrt(0.1), pw=1, lambda=1)
  }, .addProp2, addProp = 2)

  testWang2007ErrorModel("boxCox+add+pow combined 2 -> add+pow combined 2", function(f) {
    f |> model(ipre ~ add(add.sd) + pow(prop.sd, pw) + boxCox(lambda)) |>
      ini(add.sd=sqrt(0.1), prop.sd=sqrt(0.1), pw=0.5, lambda=1)
  }, c(10.886, 10.886, 10.868, 10.417, 10.662, 10.662, 10.924, 10.924,
       10.886, 10.886), addProp = 2)

  testWang2007ErrorModel("boxCox+add+pow combined 1 -> add+prop combined1", function(f) {
    f |> model(ipre ~ add(add.sd) + pow(prop.sd, pw) + boxCox(lambda)) |>
      ini(add.sd=sqrt(0.1), prop.sd=sqrt(0.1), pw=1, lambda=1)
  }, .addProp1, addProp = 1)

  .addPow1 <- c(16.231, 16.231, 16.219, 16.008, 16.093, 16.093, 16.249, 16.249,
                16.231, 16.231)

  testWang2007ErrorModel("boxCox+add+pow combined 1", function(f) {
    f |> model(ipre ~ add(add.sd) + pow(prop.sd, pw) + boxCox(lambda)) |>
      ini(add.sd=sqrt(0.1), prop.sd=sqrt(0.1), pw=0.5, lambda=1)
  }, .addPow1, addProp = 1)

  testWang2007ErrorModel("boxCox+add+pow combined 1 (override)->add+pow combined 1 (override)", function(f) {
    f |> model(ipre ~ add(add.sd) + pow(prop.sd, pw) + boxCox(lambda) + combined1()) |>
      ini(add.sd=sqrt(0.1), prop.sd=sqrt(0.1), pw=0.5, lambda=1)
  }, .addPow1, addProp = 2)

  rxode2::rxUnloadAll()
})
