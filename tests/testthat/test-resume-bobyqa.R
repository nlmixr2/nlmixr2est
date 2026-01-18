nmTest({
  test_that("resume bobyqa", {

    model <- function() {
      ini({
        TVCL <- c(0,0.01)
        TVV1 <- c(0,0.1)
        TVQ <- c(0,0.02)
        TVV2 <- c(0,0.1)
        TVR0 <- c(0,1)
        TVKON <- c(0,1)
        TVKOFF <- c(0,1)
        TVKDEG <- c(0,1)
        TVKINT <- c(0,1)
        ETA.CL ~ 0.01
        ADD.ERR.PK <- 0.1
      })
      model({
        CL <- TVCL * exp(ETA.CL)
        V1 <- TVV1
        Q <- TVQ
        V2 <- TVV2
        R0 <- TVR0
        KON <- TVKON
        KOFF <- TVKOFF
        KDEG <- TVKDEG
        KINT <- TVKINT
        KEL <-  CL/V1
        K12 <-  Q/V1
        K21 <-  Q/V2
        KD <-  KOFF/KON
        KSYN <- R0*KDEG
        A3(0) <- R0
        d/dt(A1) =  - (K12+KEL)*A1 + K21*A2 - KON*A1*A3         + KOFF*A4*V1     # nmol free drug central
        d/dt(A2) =          K12*A1 - K21*A2                                      # nmol drug peripheral
        d/dt(A3) =    KSYN - KDEG*A3          - KON*(A1/V1)*A3  + KOFF*A4        # nM free target central
        d/dt(A4) =                              KON*(A1/V1)*A3  - (KINT+KOFF)*A4 # nM complex central
        CFREE = A1/V1
        CTOT  = CFREE+A4
        CP = log(CTOT)
        CP ~ add(ADD.ERR.PK)
      })
    }
    dat <- readRDS(test_path("test-resume-bobyqa.rds"))

    d <- .nlmixr(model, dat, "focei")

    expect_true(grepl("bobyqa", d$extra))
  })
})
