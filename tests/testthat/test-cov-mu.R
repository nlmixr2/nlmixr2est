test_that("mu-ref covariates", {

  f0 <- function(){
    ini({
      lka      = log(0.9)
      lcl      = log(10.3)
      lv2      = log(1.5)
      lv3      = log(82)
      lq       = log(6)
      CONT_A_CL= 0.8
      CONT_A_V2= 1.2
      CAT_B_CL = 1.5
      prop.err <- c(0, 0.1, 1)
      eta.cl      ~ 0.1
      eta.v2      ~ 0.16
    })
    model({
      ka   =  exp(lka)
      cl   =  exp(lcl+eta.cl+CONT_A_CL*logCONT_A+CAT_B_CL*CAT_B_M)
      v2   =  exp(lv2+eta.v2+CONT_A_V2*logCONT_A)
      v3   =  exp(lv3)
      q    =  exp(lq)
      linCmt() ~ prop(prop.err)
    })
  }

  f <- nlmixr(f0)

  expect_equal(f$saemModel0,
               quote(rxModelVars({
                 ka = exp(lka)
                 cl = exp(lcl)
                 v2 = exp(lv2)
                 v3 = exp(lv3)
                 q = exp(lq)
                 rx_pred_ <- linCmt()
               })))

})
