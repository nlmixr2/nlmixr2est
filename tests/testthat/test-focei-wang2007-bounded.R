nmTest({
  skip_on_cran()
  rxode2::rxUnloadAll()
  skip_on_os("windows")

  ## logitNorm
  .logitNormAdd <- c(0.612, 0.612, 0.612, 0.612, 0.786, 0.786,
                     0.629, 0.629, 0.612, 0.612)

  testWang2007ErrorModel("logitNorm", function(f) {
    f |> model(ipre ~ logitNorm(logit.sd, 0, 12)) |>
      ini(logit.sd=sqrt(0.1))
  }, .logitNormAdd)

  .logitNormProp <- c(67.882, 67.882, 67.731, 67.765, 67.765, 67.765,
                      67.697, 67.697, 67.882, 67.882)

  testWang2007ErrorModel("logitNorm(NA)+prop", function(f) {
    f |> model(ipre ~ logitNorm(NA, 0, 12) + prop(prop.sd)) |>
      ini(prop.sd=sqrt(0.1))
  }, .logitNormProp)

  testWang2007ErrorModel("logitNorm(NA)+pow->logitNorm(NA)+prop", function(f) {
    f |> model(ipre ~ logitNorm(NA, 0, 12) + pow(prop.sd, pw)) |>
      ini(prop.sd=sqrt(0.1), pw=1)
  }, .logitNormProp)

  .logitNormPow <- c(29.055, 29.055, 29.007, 28.987, 28.989, 28.989,
                     29.012, 29.012, 29.055, 29.055)

  testWang2007ErrorModel("logitNorm(NA)+pow->logitNorm(NA)+prop", function(f) {
    f |> model(ipre ~ logitNorm(NA, 0, 12) + pow(prop.sd, pw)) |>
      ini(prop.sd=sqrt(0.1), pw=0.5)
  }, .logitNormPow)

  .logitNormAddProp1 <- c(72.699, 72.699, 72.591, 72.615, 72.615,
                          72.615, 72.554, 72.554, 72.699, 72.699)

  testWang2007ErrorModel("logitNorm+prop", function(f) {
    f |> model(ipre ~ logitNorm(logit.sd, 0, 12) + prop(prop.sd)) |>
      ini(logit.sd=sqrt(0.1), prop.sd=sqrt(0.1))
  }, .logitNormAddProp1, addProp = 1)

  testWang2007ErrorModel("logitNorm+pow->logitNorm+prop", function(f) {
    f |> model(ipre ~ logitNorm(logit.sd, 0, 12) + pow(prop.sd, pw)) |>
      ini(logit.sd=sqrt(0.1), prop.sd=sqrt(0.1), pw=1)
  }, .logitNormAddProp1, addProp = 1)

  .logitNormAddPow1 <- c(40.053, 40.053, 40.028, 40.028, 40.029,
                         40.029, 40.02, 40.02, 40.053, 40.053)

  testWang2007ErrorModel("logitNorm+prop", function(f) {
    f |> model(ipre ~ logitNorm(logit.sd, 0, 12) + pow(prop.sd, pw)) |>
      ini(logit.sd=sqrt(0.1), prop.sd=sqrt(0.1), pw=0.5)
  }, .logitNormAddPow1, addProp = 1)

  ## logitNorm + yeoJohnson
  .logitNormAddYeoJohnson <- c(9.019, 9.019, 9.019, 9.019, 9.576,
                               9.576, 9.044, 9.044, 9.019, 9.019)
  testWang2007ErrorModel("logitNorm+yeoJohnson", function(f) {
    f |> model(ipre ~ logitNorm(logit.sd, 0, 12) + yeoJohnson(lm)) |>
      ini(logit.sd=sqrt(0.1), lm=0.5)
  }, .logitNormAddYeoJohnson)

  .logitNormPropAddYeoJohnson <- c(78.136, 78.136, 77.984, 78.017,
                                   78.017, 78.017, 77.95, 77.95, 78.136, 78.136)

  testWang2007ErrorModel("logitNorm(NA)+prop+yeoJohnson", function(f) {
    f |> model(ipre ~ logitNorm(NA, 0, 12) + prop(prop.sd) + yeoJohnson(lm)) |>
      ini(prop.sd=sqrt(0.1), lm=0.5)
  }, .logitNormPropAddYeoJohnson)

  testWang2007ErrorModel("logitNorm(NA)+pow+yeoJohnson->logitNorm(NA)+prop+yeoJohnson", function(f) {
    f |> model(ipre ~ logitNorm(NA, 0, 12) + pow(prop.sd, pw) + yeoJohnson(lm)) |>
      ini(prop.sd=sqrt(0.1), lm=0.5, pw=1)
  }, .logitNormPropAddYeoJohnson)

  .logitNormPowYeoJohnson <- c(39.334, 39.334, 39.287, 39.251,
                               39.256, 39.256, 39.287, 39.287, 39.334, 39.334)

  testWang2007ErrorModel("logitNorm(NA)+pow+yeoJohnson", function(f) {
    f |> model(ipre ~ logitNorm(NA, 0, 12) + pow(prop.sd, pw) + yeoJohnson(lm)) |>
      ini(prop.sd=sqrt(0.1), lm=0.5, pw=0.5)
  }, .logitNormPowYeoJohnson)

  .logitNormAddPropAddYeoJohnson1 <- c(82.941, 82.941, 82.833,
                                       82.857, 82.857, 82.857, 82.796, 82.796, 82.941, 82.941)

  testWang2007ErrorModel("logitNorm+add+prop+yeoJohnson combined 1", function(f) {
    f |> model(ipre ~ logitNorm(logit.sd, 0, 12) + prop(prop.sd) + yeoJohnson(lm)) |>
      ini(logit.sd=sqrt(0.1), prop.sd=sqrt(0.1), lm=0.5)
  }, .logitNormAddPropAddYeoJohnson1, addProp=1)

  .logitNormAddPropAddYeoJohnson2 <- c(78.485, 78.485, 78.341,
                                       78.373, 78.373, 78.373, 78.309, 78.309, 78.485, 78.485)

  testWang2007ErrorModel("logitNorm+add+prop+yeoJohnson combined2", function(f) {
    f |> model(ipre ~ logitNorm(logit.sd, 0, 12) + prop(prop.sd) + yeoJohnson(lm)) |>
      ini(logit.sd=sqrt(0.1), prop.sd=sqrt(0.1), lm=0.5)
  }, .logitNormAddPropAddYeoJohnson2, addProp = 2)

  .logitNormAddPowAddYeoJohnson1 <- c(82.941, 82.941, 82.833, 82.857,
                                      82.857, 82.857, 82.796, 82.796, 82.941, 82.941)

  testWang2007ErrorModel("logitNorm+pow+yeoJohnson combined2", function(f) {
    f |> model(ipre ~ logitNorm(logit.sd, 0, 12) + pow(prop.sd, pw) + yeoJohnson(lm)) |>
      ini(logit.sd=sqrt(0.1), prop.sd=sqrt(0.1), lm=0.5)
  }, .logitNormAddPowAddYeoJohnson1, addProp = 1)

  .logitNormAddPowAddYeoJohnson2 <- c(78.485, 78.485, 78.341, 78.373,
                                      78.373, 78.373, 78.309, 78.309, 78.485, 78.485)

  testWang2007ErrorModel("logitNorm+pow+yeoJohnson combined2", function(f) {
    f |> model(ipre ~ logitNorm(logit.sd, 0, 12) + pow(prop.sd, pw) + yeoJohnson(lm)) |>
      ini(logit.sd=sqrt(0.1), prop.sd=sqrt(0.1), lm=0.5)
  }, .logitNormAddPowAddYeoJohnson2, addProp = 2)

  .probitNormAdd <- c(12.827, 12.827, 12.827, 12.827, 12.847, 12.847,
                      12.836, 12.836, 12.827, 12.827)

  testWang2007ErrorModel("probitNorm", function(f) {
    f |> model(ipre ~ probitNorm(logit.sd, 0, 12)) |>
      ini(logit.sd=sqrt(0.1))
  }, .probitNormAdd)

  .probitNormProp <- c(88.875, 88.875, 88.733, 88.766, 88.766,
                       88.766, 88.679, 88.679, 88.875, 88.875)

  testWang2007ErrorModel("probitNorm(NA)+prop", function(f) {
    f |> model(ipre ~ probitNorm(NA, 0, 12) + prop(prop.sd)) |>
      ini(prop.sd=sqrt(0.1))
  }, .probitNormProp)

  testWang2007ErrorModel("probitNorm(NA)+pow->probitNorm(NA)+prop", function(f) {
    f |> model(ipre ~ probitNorm(NA, 0, 12) + pow(prop.sd, pw)) |>
      ini(prop.sd=sqrt(0.1), pw=1)
  }, .probitNormProp)

  .probitNormPow <- c(48.625, 48.625, 48.579, 48.587, 48.587, 48.587,
                      48.565, 48.565, 48.625, 48.625)

  testWang2007ErrorModel("probitNorm(NA)+pow->probitNorm(NA)+prop", function(f) {
    f |> model(ipre ~ probitNorm(NA, 0, 12) + pow(prop.sd, pw)) |>
      ini(prop.sd=sqrt(0.1), pw=0.5)
  }, .probitNormPow)

  .probitNormAddProp1 <- c(93.761, 93.761, 93.661, 93.682, 93.682,
                           93.682, 93.611, 93.611, 93.761, 93.761)

  testWang2007ErrorModel("probitNorm+prop", function(f) {
    f |> model(ipre ~ probitNorm(probit.sd, 0, 12) + prop(prop.sd)) |>
      ini(probit.sd=sqrt(0.1), prop.sd=sqrt(0.1))
  }, .probitNormAddProp1, addProp = 1)

  testWang2007ErrorModel("probitNorm+pow->probitNorm+prop", function(f) {
    f |> model(ipre ~ probitNorm(probit.sd, 0, 12) + pow(prop.sd, pw)) |>
      ini(probit.sd=sqrt(0.1), prop.sd=sqrt(0.1), pw=1)
  }, .probitNormAddProp1, addProp = 1)

  .probitNormAddPow1 <- c(60.418, 60.418, 60.396, 60.401, 60.401,
                          60.401, 60.378, 60.378, 60.418, 60.418)

  testWang2007ErrorModel("probitNorm+pow, combined1", function(f) {
    f |> model(ipre ~ probitNorm(probit.sd, 0, 12) + pow(prop.sd, pw)) |>
      ini(probit.sd=sqrt(0.1), prop.sd=sqrt(0.1), pw=0.5)
  }, .probitNormAddPow1, addProp = 1)

  ## probitNorm + yeoJohnson
  .probitNormAddYeoJohnson <- c(19.69, 19.69, 19.69, 19.69, 19.729,
                                19.729, 19.699, 19.699, 19.69, 19.69)

  testWang2007ErrorModel("probitNorm+yeoJohnson", function(f) {
    f |> model(ipre ~ probitNorm(probit.sd, 0, 12) + yeoJohnson(lm)) |>
      ini(probit.sd=sqrt(0.1), lm=0.5)
  }, .probitNormAddYeoJohnson)

  .probitNormPropAddYeoJohnson <- c(96.041, 96.041, 95.899, 95.931,
                                    95.931, 95.931, 95.845, 95.845, 96.041, 96.041)

  testWang2007ErrorModel("probitNorm(NA)+prop+yeoJohnson", function(f) {
    f |> model(ipre ~ probitNorm(NA, 0, 12) + prop(prop.sd) + yeoJohnson(lm)) |>
      ini(prop.sd=sqrt(0.1), lm=0.5)
  }, .probitNormPropAddYeoJohnson)

  testWang2007ErrorModel("probitNorm(NA)+pow+yeoJohnson->probitNorm(NA)+prop+yeoJohnson", function(f) {
    f |> model(ipre ~ probitNorm(NA, 0, 12) + pow(prop.sd, pw) + yeoJohnson(lm)) |>
      ini(prop.sd=sqrt(0.1), lm=0.5, pw=1)
  }, .probitNormPropAddYeoJohnson)

  .probitNormPowYeoJohnson <- c(55.798, 55.798, 55.751, 55.758,
                                55.758, 55.758, 55.737, 55.737, 55.798, 55.798)

  testWang2007ErrorModel("probitNorm(NA)+pow+yeoJohnson", function(f) {
    f |> model(ipre ~ probitNorm(NA, 0, 12) + pow(prop.sd, pw) + yeoJohnson(lm)) |>
      ini(prop.sd=sqrt(0.1), lm=0.5, pw=0.5)
  }, .probitNormPowYeoJohnson)

  .probitNormAddPropAddYeoJohnson1 <- c(100.925, 100.925, 100.824,
                                        100.846, 100.846, 100.846, 100.774, 100.774, 100.925, 100.925)

  testWang2007ErrorModel("probitNorm+add+prop+yeoJohnson combined 1", function(f) {
    f |> model(ipre ~ probitNorm(probit.sd, 0, 12) + prop(prop.sd) + yeoJohnson(lm)) |>
      ini(probit.sd=sqrt(0.1), prop.sd=sqrt(0.1), lm=0.5)
  }, .probitNormAddPropAddYeoJohnson1, addProp=1)

  .probitNormAddPropAddYeoJohnson2 <- c(96.398, 96.398, 96.264,
                                        96.295, 96.295, 96.295, 96.213, 96.213, 96.398, 96.398)

  testWang2007ErrorModel("probitNorm+add+prop+yeoJohnson combined2", function(f) {
    f |> model(ipre ~ probitNorm(probit.sd, 0, 12) + prop(prop.sd) + yeoJohnson(lm)) |>
      ini(probit.sd=sqrt(0.1), prop.sd=sqrt(0.1), lm=0.5)
  }, .probitNormAddPropAddYeoJohnson2, addProp = 2)

  .probitNormAddPowAddYeoJohnson1 <- c(100.925, 100.925, 100.824,
                                       100.846, 100.846, 100.846, 100.774, 100.774, 100.925, 100.925)

  testWang2007ErrorModel("probitNorm+pow+yeoJohnson combined1", function(f) {
    f |> model(ipre ~ probitNorm(probit.sd, 0, 12) + pow(prop.sd, pw) + yeoJohnson(lm)) |>
      ini(probit.sd=sqrt(0.1), prop.sd=sqrt(0.1), lm=0.5)
  }, .probitNormAddPowAddYeoJohnson1, addProp = 1)

  .probitNormAddPowAddYeoJohnson2 <- c(96.398, 96.398, 96.264,
                                       96.295, 96.295, 96.295, 96.213, 96.213, 96.398, 96.398)
  testWang2007ErrorModel("probitNorm+pow+yeoJohnson combined2", function(f) {
    f |> model(ipre ~ probitNorm(probit.sd, 0, 12) + pow(prop.sd, pw) + yeoJohnson(lm)) |>
      ini(probit.sd=sqrt(0.1), prop.sd=sqrt(0.1), lm=0.5)
  }, .probitNormAddPowAddYeoJohnson2, addProp = 2)

  rxode2::rxUnloadAll()
})
