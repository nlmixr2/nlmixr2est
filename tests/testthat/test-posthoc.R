nmTest({
  test_that("test that posthoc does the correct thing with one subject", {
    skip_on_cran()

    d <- readRDS(test_path("datos_pac.rds"))

    mod.dos.cmpt <- function() {
      ini({
        # *1* Efectos fijos "theta" + Error
        # A) Parametros poblacionales y su valor poblacional
        ltcl <- log(5.42*0.001) #Ln_Cl L/Dia
        ltv1 <- log(52.4*0.001) #Ln_v1 L
        ltQ <-  log(2.26*0.001) #Ln_Q  L/Dia
        ltv2 <- log(19.6*0.001) #Ln_v2 L
        # B) Variabilidad en las observaciones
        add.err  <- 0.15 #error_aditivo
        prop.err <- 0.15 #error_proporcional
        # *2* Efectos aleatorios "omega". Variabilidad entre sujetos "omega"
        # Variabilidad de los parametros poblacionales (varianza)
        eta.cl  ~ 0.3001**2 #eta_cl
        eta.v1  ~ 0.1255**2 #eta_v1
        eta.v2  ~ 0.5165**2 #eta_v2
      })
      model({
        # First parameters are defined in terms of the initial estimates
        # parameter names.
        cl <- exp(ltcl+LWT -0.313*LW65 -0.855*LALB41+ ATI*log(1.292)+ IMM*log(0.863)  + eta.cl)
        v1 <- exp(ltv1+LWT -0.233*LW65 +eta.v1)
        Q  <- exp(ltQ+LWT)
        v2 <- exp(ltv2 +LWT-0.588*LW65  +eta.v2)
        # After the differential equations are defined
        k <- cl/v1
        k12 <- Q/v1
        k21 <- Q/v2
        d/dt(central)  = -(k+k12)*central+k21*peripheral # masa (mg) de IFX en compart. central
        d/dt(peripheral) = k12*central-k21*peripheral    #Free Drug second compartment amount
        # And the concentration is then calculated
        conc=central/v1
        conc ~ add(add.err) + prop(prop.err)
      })
    }
    # Calculations required by the model of interest
    funcion_transf_logaritmica<-function(nombre_dataframe){
      nombre_dataframe$LWT    <-  log(nombre_dataframe$WT)
      nombre_dataframe$LW65   <-  log(nombre_dataframe$WT/65)
      nombre_dataframe$LALB41 <-  log(nombre_dataframe$ALB/4.1)
      nombre_dataframe <- data.frame(nombre_dataframe)
    }

    d <- funcion_transf_logaritmica(d)

    f <- .nlmixr(mod.dos.cmpt, d, "posthoc", control=list(calcTables=FALSE))

    f2 <- suppressMessages(addTable(f))

    expect_false(all(f2$eta.cl == 0.0))
    expect_false(all(f2$eta.v1 == 0.0))
    expect_false(all(f2$eta.v2 == 0.0))

  })
})
