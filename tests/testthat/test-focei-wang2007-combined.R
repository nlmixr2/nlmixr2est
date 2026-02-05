nmTest({
  skip_on_cran()
  rxode2::rxUnloadAll()
  skip_on_os("windows")

  # Expected values for combined error model tests
  .addVals <- c(-2.059, -2.059, -2.059, -2.059, 0.026, 0.026, -1.997,
                -1.997, -2.059, -2.059)
  .addModVals <- c(3.238, 3.238, 3.207, 2.438, 3.311, 3.311, 3.298,
                   3.298, 3.238, 3.238)
  .propVals <- c(39.458, 39.458, 39.275, 39.207, 39.213, 39.213,
                 39.365, 39.365, 39.458, 39.458)
  .propModVals <- c(63.353, 63.353, 63.001, 63.063, 63.063, 63.063,
                    63.063, 63.063, 63.353, 63.353)
  .propFVals <- c(6.496, 6.496, 6.488, 6.275, 9.262, 9.262, 6.545,
                  6.545, 6.496, 6.496)
  .propFModVals <- c(19.177, 19.177, 19.07, 18.202, 18.333, 18.333,
                     19.158, 19.158, 19.177, 19.177)
  .addModPropMod2 <- c(106.308, 106.308, 106.013, 106.079, 106.079,
                       106.079, 105.948, 105.948, 106.308, 106.308)
  .addProp2 <- c(39.735, 39.735, 39.562, 39.499, 39.505, 39.505,
                 39.647, 39.647, 39.735, 39.735)
  .addModPropMod1 <- c(106.308, 106.308, 106.013, 106.079, 106.079,
                       106.079, 105.948, 105.948, 106.308, 106.308)
  .addProp1 <- c(43.554, 43.554, 43.416, 43.394, 43.398, 43.398,
                 43.469, 43.469, 43.554, 43.554)
  .addModPropFModVals2 <- c(54.317, 54.317, 54.14, 54.165, 54.166,
                            54.166, 54.148, 54.148, 54.317, 54.317)
  .addPropFVals2 <- c(-2.321, -2.321, -2.322, -2.454, -0.65, -0.65,
                      -2.247, -2.247, -2.321, -2.321)
  .addPow1 <- c(16.231, 16.231, 16.219, 16.008, 16.093, 16.093, 16.249, 16.249,
                16.231, 16.231)

  ################################################################################
  # Add+Proportional tests (combined 2)
  ################################################################################
  testWang2007ErrorModel("add+prop, combined 2->add", function(f) {
    f |> model(ipre ~ add(add.sd) + prop(prop.sd)) |> ini(add.sd=sqrt(0.1), prop.sd=0)
  }, .addVals, addProp = 2)

  testWang2007ErrorModel("addMod+prop, combined 2->addMod", function(f) {
    f |> model(ipre ~ add(f2) + prop(prop.sd)) |> ini(prop.sd=0)
  }, .addModVals, addProp = 2)

  testWang2007ErrorModel("add+prop, combined 2->prop", function(f) {
    f |> model(ipre ~ add(add.sd) + prop(prop.sd)) |> ini(add.sd=0, prop.sd=sqrt(0.1))
  }, .propVals, addProp = 2)

  testWang2007ErrorModel("add+propMod, combined 2->propMod", function(f) {
    f |> model(ipre ~ add(add.sd) + prop(f2)) |> ini(add.sd=0)
  }, .propModVals, addProp = 2)

  testWang2007ErrorModel("addMod+propMod, combined 2->propMod", function(f) {
    f |> model(ipre ~ add(f2) + prop(f3))
  },  .addModPropMod2, addProp = 2)

  testWang2007ErrorModel("add+prop, combined 2", function(f) {
    f |> model(ipre ~ add(add.sd) + prop(prop.sd)) |> ini(add.sd=sqrt(0.1), prop.sd=sqrt(0.1))
  }, .addProp2, addProp = 2)

  testWang2007ErrorModel("add+prop, combined 2", function(f) {
    f |> model(ipre ~ add(add.sd) + prop(prop.sd) + combined2()) |> ini(add.sd=sqrt(0.1), prop.sd=sqrt(0.1))
  }, .addProp2, addProp = 1)

  ## propT
  testWang2007ErrorModel("add+propT, combined 2->add", function(f) {
    f |> model(ipre ~ add(add.sd) + propT(prop.sd)) |> ini(add.sd=sqrt(0.1), prop.sd=0)
  }, .addVals, addProp = 2)

  testWang2007ErrorModel("addMod+propT, combined 2->addMod", function(f) {
    f |> model(ipre ~ add(f2) + propT(prop.sd)) |> ini(prop.sd=0)
  }, .addModVals, addProp = 2)

  testWang2007ErrorModel("add+propT, combined 2->prop", function(f) {
    f |> model(ipre ~ add(add.sd) + propT(prop.sd)) |> ini(add.sd=0, prop.sd=sqrt(0.1))
  }, .propVals, addProp = 2)

  testWang2007ErrorModel("add+propTMod, combined 2->propMod", function(f) {
    f |> model(ipre ~ add(add.sd) + propT(f2)) |> ini(add.sd=0)
  }, .propModVals, addProp = 2)

  testWang2007ErrorModel("addMod+propTMod, combined 2->propMod", function(f) {
    f |> model(ipre ~ add(f2) + propT(f3))
  },  .addModPropMod2, addProp = 2)

  testWang2007ErrorModel("add+propT, combined 2", function(f) {
    f |> model(ipre ~ add(add.sd) + propT(prop.sd)) |> ini(add.sd=sqrt(0.1), prop.sd=sqrt(0.1))
  }, .addProp2, addProp = 2)

  testWang2007ErrorModel("add+propT, combined 2 (specified)", function(f) {
    f |> model(ipre ~ add(add.sd) + propT(prop.sd) + combined2()) |> ini(add.sd=sqrt(0.1), prop.sd=sqrt(0.1))
  }, .addProp2, addProp = 1)

  # propF
  testWang2007ErrorModel("add+propF, combined 2->add", function(f) {
    f |> model(ipre ~ add(add.sd) + propF(prop.sd, f2)) |> ini(add.sd=sqrt(0.1), prop.sd=0)
  }, .addVals, addProp = 2)

  testWang2007ErrorModel("addMod+propF, combined 2->addMod", function(f) {
    f |> model(ipre ~ add(f2) + propF(prop.sd, f2)) |> ini(prop.sd=0)
  }, .addModVals, addProp = 2)

  testWang2007ErrorModel("add+propF, combined 2->prop", function(f) {
    f |> model(ipre ~ add(add.sd) + propF(prop.sd, f2)) |> ini(add.sd=0, prop.sd=sqrt(0.1))
  }, .propFVals, addProp = 2)

  testWang2007ErrorModel("add+propFMod, combined 2->propMod", function(f) {
    f |> model(ipre ~ add(add.sd) + propF(lipre, f2)) |> ini(add.sd=0)
  }, .propFModVals, addProp = 2)

  testWang2007ErrorModel("addMod+propFMod, combined 2->propMod", function(f) {
    f |> model(ipre ~ add(f2) + propF(lipre, f3))
  }, .addModPropFModVals2, addProp = 2)

  testWang2007ErrorModel("add+propF, combined 2", function(f) {
    f |> model(ipre ~ add(add.sd) + propF(prop.sd, f2)) |> ini(add.sd=sqrt(0.1), prop.sd=sqrt(0.1))
  },  .addPropFVals2, addProp = 2)

  testWang2007ErrorModel("add+propF, combined 2 (specified)", function(f) {
    f |> model(ipre ~ add(add.sd) + propF(prop.sd, f2) + combined2()) |> ini(add.sd=sqrt(0.1), prop.sd=sqrt(0.1))
  }, .addPropFVals2, addProp = 1)

  ################################################################################
  # Add+Proportional tests (combined 1)
  ################################################################################
  testWang2007ErrorModel("add+prop, combined 1->add", function(f) {
    f |> model(ipre ~ add(add.sd) + prop(prop.sd)) |> ini(add.sd=sqrt(0.1), prop.sd=0)
  }, .addVals, addProp = 1)

  testWang2007ErrorModel("addMod+prop, combined 1->addMod", function(f) {
    f |> model(ipre ~ add(f2) + prop(prop.sd)) |> ini(prop.sd=0)
  }, .addModVals, addProp = 1)

  testWang2007ErrorModel("add+prop, combined 1->prop", function(f) {
    f |> model(ipre ~ add(add.sd) + prop(prop.sd)) |> ini(add.sd=0, prop.sd=sqrt(0.1))
  }, .propVals, addProp = 1)

  testWang2007ErrorModel("add+propMod, combined 1->propMod", function(f) {
    f |> model(ipre ~ add(add.sd) + prop(f2)) |> ini(add.sd=0)
  }, .propModVals, addProp = 1)

  testWang2007ErrorModel("addMod+propMod, combined 1->propMod", function(f) {
    f |> model(ipre ~ add(f2) + prop(f3))
  },  .addModPropMod1, addProp = 2)

  testWang2007ErrorModel("add+prop, combined 2", function(f) {
    f |> model(ipre ~ add(add.sd) + prop(prop.sd)) |> ini(add.sd=sqrt(0.1), prop.sd=sqrt(0.1))
  }, .addProp2, addProp = 2)

  testWang2007ErrorModel("add+prop, combined 2", function(f) {
    f |> model(ipre ~ add(add.sd) + prop(prop.sd) + combined2()) |> ini(add.sd=sqrt(0.1), prop.sd=sqrt(0.1))
  }, .addProp2, addProp = 1)

  ## propT
  testWang2007ErrorModel("add+propT, combined 2->add", function(f) {
    f |> model(ipre ~ add(add.sd) + propT(prop.sd)) |> ini(add.sd=sqrt(0.1), prop.sd=0)
  }, .addVals, addProp = 2)

  testWang2007ErrorModel("addMod+propT, combined 2->addMod", function(f) {
    f |> model(ipre ~ add(f2) + propT(prop.sd)) |> ini(prop.sd=0)
  }, .addModVals, addProp = 2)

  testWang2007ErrorModel("add+propT, combined 2->prop", function(f) {
    f |> model(ipre ~ add(add.sd) + propT(prop.sd)) |> ini(add.sd=0, prop.sd=sqrt(0.1))
  }, .propVals, addProp = 2)

  testWang2007ErrorModel("add+propTMod, combined 2->propMod", function(f) {
    f |> model(ipre ~ add(add.sd) + propT(f2)) |> ini(add.sd=0)
  }, .propModVals, addProp = 2)

  testWang2007ErrorModel("addMod+propTMod, combined 2->propMod", function(f) {
    f |> model(ipre ~ add(f2) + propT(f3))
  },  .addModPropMod2, addProp = 2)

  testWang2007ErrorModel("add+propT, combined 2", function(f) {
    f |> model(ipre ~ add(add.sd) + propT(prop.sd)) |> ini(add.sd=sqrt(0.1), prop.sd=sqrt(0.1))
  }, .addProp2, addProp = 2)

  testWang2007ErrorModel("add+propT, combined 2 (specified)", function(f) {
    f |> model(ipre ~ add(add.sd) + propT(prop.sd) + combined2()) |> ini(add.sd=sqrt(0.1), prop.sd=sqrt(0.1))
  }, .addProp2, addProp = 1)

  # propF
  testWang2007ErrorModel("add+propF, combined 2->add", function(f) {
    f |> model(ipre ~ add(add.sd) + propF(prop.sd, f2)) |> ini(add.sd=sqrt(0.1), prop.sd=0)
  }, .addVals, addProp = 2)

  testWang2007ErrorModel("addMod+propF, combined 2->addMod", function(f) {
    f |> model(ipre ~ add(f2) + propF(prop.sd, f2)) |> ini(prop.sd=0)
  }, .addModVals, addProp = 2)

  testWang2007ErrorModel("add+propF, combined 2->prop", function(f) {
    f |> model(ipre ~ add(add.sd) + propF(prop.sd, f2)) |> ini(add.sd=0, prop.sd=sqrt(0.1))
  }, .propFVals, addProp = 2)

  testWang2007ErrorModel("add+propFMod, combined 2->propMod", function(f) {
    f |> model(ipre ~ add(add.sd) + propF(lipre, f2)) |> ini(add.sd=0)
  }, .propFModVals, addProp = 2)

  testWang2007ErrorModel("addMod+propFMod, combined 2->propMod", function(f) {
    f |> model(ipre ~ add(f2) + propF(lipre, f3))
  }, .addModPropFModVals2, addProp = 2)

  testWang2007ErrorModel("add+propF, combined 2", function(f) {
    f |> model(ipre ~ add(add.sd) + propF(prop.sd, f2)) |> ini(add.sd=sqrt(0.1), prop.sd=sqrt(0.1))
  },  .addPropFVals2, addProp = 2)

  testWang2007ErrorModel("add+propF, combined 2 (specified)", function(f) {
    f |> model(ipre ~ add(add.sd) + propF(prop.sd, f2) + combined2()) |> ini(add.sd=sqrt(0.1), prop.sd=sqrt(0.1))
  }, .addPropFVals2, addProp = 1)

  #################################################################################################
  testWang2007ErrorModel("add+prop, combined 1", function(f) {
    f |> model(ipre ~ add(add.sd) + prop(prop.sd)) |> ini(add.sd=sqrt(0.1), prop.sd=sqrt(0.1))
  }, .addProp1, addProp = 1)

  testWang2007ErrorModel("add+prop, combined 1->add", function(f) {
    f |> model(ipre ~ add(add.sd) + prop(prop.sd)) |> ini(add.sd=sqrt(0.1), prop.sd=0)
  }, .addVals, addProp = 1)

  testWang2007ErrorModel("add+prop, combined 1->prop", function(f) {
    f |> model(ipre ~ add(add.sd) + prop(prop.sd)) |> ini(add.sd=0, prop.sd=sqrt(0.1))
  }, .propVals, addProp = 1)

  testWang2007ErrorModel("add+prop, combined 1 (override)", function(f) {
    f |> model(ipre ~ add(add.sd) + prop(prop.sd) + combined1()) |> ini(add.sd=sqrt(0.1), prop.sd=sqrt(0.1))
  }, .addProp1, addProp = 2)

  testWang2007ErrorModel("add+prop, combined 2 (override)", function(f) {
    f |> model(ipre ~ add(add.sd) + prop(prop.sd) + combined2()) |> ini(add.sd=sqrt(0.1), prop.sd=sqrt(0.1))
  }, .addProp2, addProp = 1)

  ################################################################################
  # Add+Pow tests (combined 2)
  ################################################################################

  testWang2007ErrorModel("add+pow combined 2 -> add+prop combined2", function(f) {
    f |> model(ipre ~ add(add.sd) + pow(prop.sd, pw)) |> ini(add.sd=sqrt(0.1), prop.sd=sqrt(0.1), pw=1)
  }, .addProp2, addProp = 2)

  testWang2007ErrorModel("add+pow combined 2", function(f) {
    f |> model(ipre ~ add(add.sd) + pow(prop.sd, pw)) |> ini(add.sd=sqrt(0.1), prop.sd=sqrt(0.1), pw=0.5)
  }, c(10.886, 10.886, 10.868, 10.417, 10.662, 10.662, 10.924, 10.924,
       10.886, 10.886), addProp = 2)

  testWang2007ErrorModel("add+pow combined 1->add+prop combined1", function(f) {
    f |> model(ipre ~ add(add.sd) + pow(prop.sd, pw)) |> ini(add.sd=sqrt(0.1), prop.sd=sqrt(0.1), pw=1)
  }, .addProp1, addProp = 1)

  testWang2007ErrorModel("add+pow combined 1", function(f) {
    f |> model(ipre ~ add(add.sd) + pow(prop.sd, pw)) |> ini(add.sd=sqrt(0.1), prop.sd=sqrt(0.1), pw=0.5)
  }, .addPow1, addProp = 1)

  testWang2007ErrorModel("add+pow combined 1 (override)", function(f) {
    f |> model(ipre ~ add(add.sd) + pow(prop.sd, pw) + combined1()) |> ini(add.sd=sqrt(0.1), prop.sd=sqrt(0.1), pw=0.5)
  }, .addPow1, addProp = 2)

  rxode2::rxUnloadAll()
})
