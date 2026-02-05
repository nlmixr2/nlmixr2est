nmTest({
  
  # Setup data for NONMEM comparison tests
  dat <- getWang2007BaseData()
  dat2 <- getWang2007DoseData()
  
  # Get base models for NONMEM comparison
  f <- getWang2007BaseModel()
  fo <- getWang2007OdeModel()
  
  fit.prop <- .nlmixr(f, dat, "focei", foceiControl(maxOuterIterations = 0, covMethod = ""))
  fit.prop2 <- .nlmixr(fo, dat2, "focei", foceiControl(maxOuterIterations = 0, covMethod = ""))
  
  out.focei.prop <- readRDS(test_path("out.focei.prop.rds"))
  
  test_that("Matches NONMEM objective proportional function; (Based on Wang2007)", {
    expect_equal(fit.prop$objective, 39.458, tolerance=1e-3) # Matches Table 2 Prop FOCEI for NONMEM
    expect_equal(fit.prop$eta.ke, out.focei.prop$ETA1, tolerance=1e-3) # match NONMEM output
    ## Individual properties
    expect_equal(fit.prop$IPRED, out.focei.prop$IPRE, tolerance=1e-3)
    expect_equal(fit.prop$IRES, out.focei.prop$IRES, tolerance=1e-3)
    expect_equal(fit.prop$IWRES, out.focei.prop$IWRES, tolerance=1e-3)
    ## WRES variants
    expect_equal(fit.prop$PRED, out.focei.prop$NPRED, tolerance=1e-3) # matches output of PRED from NONMEM
    expect_equal(fit.prop$PRED, out.focei.prop$PRED, tolerance=1e-3) # matches output of PRED from NONMEM
    expect_equal(fit.prop$RES, out.focei.prop$RES, tolerance=1e-3) # match NONMEM output
    expect_equal(fit.prop$RES, out.focei.prop$NRES, tolerance=1e-3) # match NONMEM output
    ## FOI equivalents
    expect_equal(fit.prop$PRED, out.focei.prop$PREDI, tolerance=1e-3) # matches output of PRED from NONMEM
    ## CWRES variants
    expect_equal(fit.prop$CRES, out.focei.prop$CRES, tolerance=1e-3) # match NONMEM output
    expect_equal(fit.prop$CPRED, out.focei.prop$CPRED, tolerance=1e-3) # match NONMEM output
    expect_equal(fit.prop$CWRES, out.focei.prop$CWRES, tolerance=1e-3) # match NONMEM output
    ## Note that E[x] for CPRED and CPREDI are equal
    expect_equal(fit.prop$CRES, out.focei.prop$CRESI, tolerance=1e-3) # match NONMEM output
    expect_equal(fit.prop$CPRED, out.focei.prop$CPREDI, tolerance=1e-3) # match NONMEM output
  })
  
  
  test_that("Matches NONMEM objective proportional function; (Based on Wang2007; unoptimized)", {
    # Check unoptimized expression
    expect_equal(fit.prop2$objective, 39.458, tolerance=1e-3) # Matches Table 2 Prop FOCEI for NONMEM
    expect_equal(fit.prop2$eta.ke, out.focei.prop$ETA1, tolerance=1e-3) # match NONMEM output
    ## Individual properties
    expect_equal(fit.prop2$IPRED, out.focei.prop$IPRE, tolerance=1e-3)
    expect_equal(fit.prop2$IRES, out.focei.prop$IRES, tolerance=1e-3)
    expect_equal(fit.prop2$IWRES, out.focei.prop$IWRES, tolerance=1e-3)
    ## WRES variants
    expect_equal(fit.prop2$PRED, out.focei.prop$NPRED, tolerance=1e-3) # matches output of PRED from NONMEM
    expect_equal(fit.prop2$PRED, out.focei.prop$PRED, tolerance=1e-3) # matches output of PRED from NONMEM
    expect_equal(fit.prop2$RES, out.focei.prop$RES, tolerance=1e-3) # match NONMEM output
    expect_equal(fit.prop2$RES, out.focei.prop$NRES, tolerance=1e-3) # match NONMEM output
    ## FOI equivalents
    expect_equal(fit.prop2$PRED, out.focei.prop$PREDI, tolerance=1e-3) # matches output of PRED from NONMEM
    ## CWRES variants
    expect_equal(fit.prop2$CRES, out.focei.prop$CRES, tolerance=1e-3) # match NONMEM output
    expect_equal(fit.prop2$CPRED, out.focei.prop$CPRED, tolerance=1e-3) # match NONMEM output
    expect_equal(fit.prop2$CWRES, out.focei.prop$CWRES, tolerance=1e-3) # match NONMEM output
    ## Note that E[x] for CPRED and CPREDI are equal
    expect_equal(fit.prop2$CRES, out.focei.prop$CRESI, tolerance=1e-3) # match NONMEM output
    expect_equal(fit.prop2$CPRED, out.focei.prop$CPREDI, tolerance=1e-3) # match NONMEM output
  })
  
  # Expected value vectors for basic tests
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
  
  ################################################################################
  # Proportional tests
  ################################################################################
  testWang2007ErrorModel("prop", function(f) {
    f |> model(ipre ~ prop(prop.sd)) |> ini(prop.sd=sqrt(0.1))
  }, .propVals)
  
  testWang2007ErrorModel("propMod", function(f) {
    f |> model(ipre ~ prop(f2))
  }, .propModVals)
  
  # In this case propT = prop
  testWang2007ErrorModel("propT", function(f) {
    f |> model(ipre ~ propT(prop.sd)) |> ini(prop.sd=sqrt(0.1))
  }, .propVals)
  
  testWang2007ErrorModel("propTMod", function(f) {
    f |> model(ipre ~ propT(f2))
  }, .propModVals)
  
  testWang2007ErrorModel("propF", function(f) {
    f |> model(ipre ~ propF(prop.sd, f2)) |> ini(prop.sd=sqrt(0.1))
  }, .propFVals)
  
  testWang2007ErrorModel("propFMod", function(f) {
    f |> model(ipre ~ propF(lipre, f2))
  }, .propFModVals)
  
  ################################################################################
  # Additive Model Tests
  ################################################################################
  testWang2007ErrorModel("add", function(f) {
    f |> model(ipre ~ add(add.sd)) |> ini(add.sd=sqrt(0.1))
  }, .addVals)
  
  testWang2007ErrorModel("addMod", function(f) {
    f |> model(ipre ~ add(f2))
  }, .addModVals)
  
  skip_on_cran()
  rxode2::rxUnloadAll() # don't do too much on windows because of dll overloading
  skip_on_os("windows")
})
