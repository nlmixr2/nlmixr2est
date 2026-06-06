nmTest({
  skip_on_cran()
  rxode2::rxUnloadAll()
  skip_on_os("windows")

  # Now yeoJohnson
  .yeoJohnsonAdd <- c(11.279, 11.279, 11.279, 11.279, 11.738, 11.738,
                      11.32, 11.32, 11.279, 11.279)

  testWang2007ErrorModel("add+yeoJohnson", function(f) {
    f |> model(ipre ~ add(add.sd) + yeoJohnson(lm)) |> ini(add.sd=sqrt(0.1), lm=0.5)
  }, .yeoJohnsonAdd)

  .yeoJohnsonProp <- c(80.257, 80.257, 80.102, 80.136, 80.136,
                       80.136, 80.076, 80.076, 80.257, 80.257)
  testWang2007ErrorModel("prop+yeoJohnson", function(f) {
    f |> model(ipre ~ prop(prop.sd) + yeoJohnson(lm)) |> ini(prop.sd=sqrt(0.1), lm=0.5)
  }, .yeoJohnsonProp)

  testWang2007ErrorModel("pow+yeoJohnson", function(f) {
    f |> model(ipre ~ pow(prop.sd, pw) + yeoJohnson(lm)) |> ini(prop.sd=sqrt(0.1), lm=0.5, pw=0.5)
  }, c(41.644, 41.644, 41.598, 41.552, 41.559, 41.559, 41.607, 41.607,
       41.644, 41.644))

  .yeoJohnsonAddProp1 <- c(85.048, 85.048, 84.937, 84.961, 84.962,
                           84.962, 84.905, 84.905, 85.048, 85.048)

  testWang2007ErrorModel("add+prop+yeoJohnson", function(f) {
    f |> model(ipre ~ add(add.sd) + prop(prop.sd) + yeoJohnson(lm)) |>
      ini(add.sd=sqrt(0.1), prop.sd=sqrt(0.1), lm=0.5)
  }, .yeoJohnsonAddProp1, addProp = 1)

  .yeoJohnsonAddProp2 <- c(80.605, 80.605, 80.458, 80.49, 80.491,
                           80.491, 80.433, 80.433, 80.605, 80.605)

  testWang2007ErrorModel("add+prop+yeoJohnson", function(f) {
    f |> model(ipre ~ add(add.sd) + prop(prop.sd) + yeoJohnson(lm)) |>
      ini(add.sd=sqrt(0.1), prop.sd=sqrt(0.1), lm=0.5)
  }, .yeoJohnsonAddProp2, addProp = 2)

  .yeoJohnsonAddPow1 <- c(52.481, 52.481, 52.456, 52.451, 52.452,
                          52.452, 52.451, 52.451, 52.481, 52.481)

  testWang2007ErrorModel("add+pow+yeoJohnson", function(f) {
    f |> model(ipre ~ add(add.sd) + pow(prop.sd, pw) + yeoJohnson(lm)) |>
      ini(add.sd=sqrt(0.1), prop.sd=sqrt(0.1), pw=0.5, lm=0.5)
  }, .yeoJohnsonAddPow1, addProp = 1)

  .yeoJohnsonAddPow2 <- c(43.704, 43.704, 43.669, 43.641, 43.645,
                          43.645, 43.674, 43.674, 43.704, 43.704)

  testWang2007ErrorModel("add+pow+yeoJohnson", function(f) {
    f |> model(ipre ~ add(add.sd) + pow(prop.sd, pw) + yeoJohnson(lm)) |>
      ini(add.sd=sqrt(0.1), prop.sd=sqrt(0.1), pw=0.5, lm=0.5)
  }, .yeoJohnsonAddPow2, addProp = 2)

  rxode2::rxUnloadAll()
})
