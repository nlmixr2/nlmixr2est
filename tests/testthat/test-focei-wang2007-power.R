nmTest({
  
  skip_on_cran()
  rxode2::rxUnloadAll()
  skip_on_os("windows")
  
  # Expected value vectors for power tests
  .propVals <- c(39.458, 39.458, 39.275, 39.207, 39.213, 39.213,
                 39.365, 39.365, 39.458, 39.458)
  .propModVals <- c(63.353, 63.353, 63.001, 63.063, 63.063, 63.063,
                    63.063, 63.063, 63.353, 63.353)
  .propFVals <- c(6.496, 6.496, 6.488, 6.275, 9.262, 9.262, 6.545,
                  6.545, 6.496, 6.496)
  .powVals <- c(9.966, 9.966, 9.948, 9.331, 9.651, 9.651, 10.007,
                10.007, 9.966, 9.966)
  .powFModVals <- c(0.776, 0.776, 0.772, 0.58, 3.2, 3.2, 0.828,
                    0.828, 0.776, 0.776)
  .powF1Vals <- c(27.301, 27.301, 27.147, 26.854, 26.888, 26.888,
                  27.259, 27.259, 27.301, 27.301)
  .powF2Vals <- c(17.831, 17.831, 17.785, 16.617, 16.852, 16.852,
                  17.871, 17.871, 17.831, 17.831)
  .powF3Vals <- c(79.733, 79.733, 79.371, 79.448, 79.448, 79.448,
                  79.389, 79.389, 79.733, 79.733)
  .powFMod1Vals <- c(24.877, 24.877, 24.773, 24.518, 24.554, 24.554,
                     24.842, 24.842, 24.877, 24.877)
  .powFMod2Vals <- c(49.312, 49.312, 49.311, 49.296, 52.603, 52.603,
                     49.34, 49.34, 49.312, 49.312)
  .powFMod3Vals <- c(10.848, 10.848, 10.784, 9.446, 9.96, 9.96,
                     10.915, 10.915, 10.848, 10.848)
  
  ################################################################################
  # Power Model Tests
  ################################################################################
  testWang2007ErrorModel("pow1=prop", function(f) {
    f |> model(ipre ~ pow(pow.sd, pw)) |> ini(pow.sd=sqrt(0.1), pw=1)
  }, .propVals)
  
  testWang2007ErrorModel("powT1=propT=prop", function(f) {
    f |> model(ipre ~ powT(pow.sd, pw)) |> ini(pow.sd=sqrt(0.1), pw=1)
  }, .propVals)
  
  testWang2007ErrorModel("powF1=propF", function(f) {
    f |> model(ipre ~ powF(pow.sd, pw, f2)) |> ini(pow.sd=sqrt(0.1), pw=1)
  }, .propFVals)
  
  testWang2007ErrorModel("pow1Mod=propMod", function(f) {
    f |> model(ipre ~ pow(f2, pw)) |> ini(pw=1)
  }, .propModVals)
  
  testWang2007ErrorModel("pow", function(f) {
    f |> model(ipre ~ pow(pow.sd, pw)) |> ini(pow.sd=sqrt(0.1), pw=0.5)
  }, .powVals)
  
  testWang2007ErrorModel("powF1", function(f) {
    f |> model(ipre ~ pow(f2, pw)) |> ini(pw=0.5)
  }, .powF1Vals)
  
  testWang2007ErrorModel("powF2", function(f) {
    f |> model(ipre ~ pow(pow.sd, f2)) |> ini(pow.sd=sqrt(0.1))
  }, .powF2Vals)
  
  testWang2007ErrorModel("powF3", function(f) {
    f |> model(ipre ~ pow(lipre, f2))
  }, .powF3Vals)
  
  testWang2007ErrorModel("powT", function(f) {
    f |> model(ipre ~ powT(pow.sd, pw)) |> ini(pow.sd=sqrt(0.1), pw=0.5)
  }, .powVals)
  
  testWang2007ErrorModel("powTF1", function(f) {
    f |> model(ipre ~ powT(f2, pw)) |> ini(pw=0.5)
  }, .powF1Vals)
  
  testWang2007ErrorModel("powFT2", function(f) {
    f |> model(ipre ~ powT(pow.sd, f2)) |> ini(pow.sd=sqrt(0.1))
  }, .powF2Vals)
  
  testWang2007ErrorModel("powTF3", function(f) {
    f |> model(ipre ~ powT(lipre, f2))
  }, .powF3Vals)
  
  testWang2007ErrorModel("powFMod", function(f) {
    f |> model(ipre ~ powF(pow.sd, pw, f2)) |> ini(pow.sd=sqrt(0.1), pw=0.5)
  }, .powFModVals)
  
  testWang2007ErrorModel("powFMod1", function(f) {
    f |> model(ipre ~ powF(lipre, pw, f2)) |> ini(pw=0.5)
  }, .powFMod1Vals)
  
  testWang2007ErrorModel("powFMod2", function(f) {
    f |> model(ipre ~ powF(pow.sd, lipre, f2)) |> ini(pow.sd=sqrt(0.1))
  }, .powFMod2Vals)
  
  testWang2007ErrorModel("powFMod3", function(f) {
    f |> model(ipre ~ powF(lipre, f3, f2))
  }, .powFMod3Vals)
  
  rxode2::rxUnloadAll()
})
