nmTest({

  test_that("test augPred", {

    PKdata <- nlmixr2data::warfarin %>%
      dplyr::filter(dvid == "cp") %>%
      dplyr::select(-dvid) %>%
      dplyr::mutate(sex = ifelse(sex == "male", 1, 0))

    One.comp.KA.solved <- function() {
      ini({
        # Where initial conditions/variables are specified
        lka  <- log(1.15)  #log ka (1/h)
        lcl  <- log(0.135) #log Cl (L/h)
        lv   <- log(8)     #log V (L)
        prop.err <- 0.15   #proportional error (SD/mean)
        add.err  <- 0.6    #additive error (mg/L)
        eta.ka ~ 0.5   #IIV ka
        eta.cl ~ 0.1   #IIV cl
        eta.v  ~ 0.1   #IIV v
      })
      model({
        # Where the model is specified
        cl <- exp(lcl + eta.cl)
        v  <- exp(lv + eta.v)
        ka <- exp(lka + eta.ka)
        linCmt() ~ prop(prop.err) + add(add.err)
      })
    }


    fitOne.comp.KA.solved_S <-
      nlmixr(
        One.comp.KA.solved,    #the model definition
        PKdata,                #the data set
        est = "saem",          #the estimation algorithm (SAEM)
        #the SAEM minimisation options:
        saemControl(nBurn = 200, #200 SAEM burn-in iterations (the default)
                    nEm   = 300, #300 EM iterations (the default)
                    print = 50), #print every 50th iteration
        #only print every 50th estimation step (default=1 which gives endless output)
        tableControl(cwres = TRUE,npde=TRUE) #calculates NONMEM-style conditional weighted residuals and npde for diagnostics
      )

      expect_error(augPred(fitOne.comp.KA.solved_S), NA)

      # Data' purpose illustrates the error and my data set
      df <- tibble::tibble(
        ID = c(rep(1, 6), rep(2, 6)),
        TIME = c(0.00, 12.11, 18.41, 23.89, 36.00, 43.51, 0.00, 12.00, 20.00, 24.00, 36.80, 45.00),
        AMT = c(1000, 1000, NA, 1000, 1000, NA, 1000, 2000, NA, 1000, 1000, NA),
        DUR = c(2.5, 2.5, NA, 2.5, 2.5, NA, 2.5, 2.5, NA, 2.5, 2.5, NA),
        DV = c(NA, NA, 3.0, NA, NA, 9.6, NA, NA, 7.0, NA, NA, 2.8),
        WT = c(rep(55, 6), rep(48, 6))
      ) %>%
        dplyr::mutate(EVID = ifelse(is.na(DV), 1, 0))

      fun <- function() {
        ini({
          tvCl <- c(0, 4, Inf)
          tvVc <- c(0, 48, Inf)

          eta.Vc ~ 0.62
          prop.sd <- 0.051529

        })
        model({
          Cl <- tvCl
          Vc <- tvVc*(WT/70)*exp(eta.Vc)

          # dynamical system
          linCmt() ~ prop(prop.sd)
        })
      }

      fit <- nlmixr2(fun, df, list(print=0), est="posthoc")

      expect_error(augPred(fit), NA)

  })

  test_that("test augPred with xgxr dataset", {

    dat <- xgxr::case1_pkpd %>%
      dplyr::rename(DV=LIDV) %>%
      dplyr::filter(CMT %in% 1:2) %>%
      dplyr::filter(TRTACT != "Placebo")

      doses <- unique(dat$DOSE)
      nid <- 3 # 7 ids per dose group
      dat2 <- do.call("rbind",
                      lapply(doses, function(x) {
                        ids <- dat %>%
                          dplyr::filter(DOSE == x) %>%
                          dplyr::summarize(ids=unique(ID)) %>%
                          dplyr::pull()
                        ids <- ids[seq(1, nid)]
                        dat %>%
                          dplyr::filter(ID %in% ids)
                      }))

      cmt2 <- function(){
        ini({
          lka <- log(0.1) # log Ka
          lv <- log(10) # Log Vc
          lcl <- log(4) # Log Cl
          lq <- log(10) # log Q
          lvp <- log(20) # Log Vp

          eta.ka ~ 0.01
          eta.v ~ 0.1
          eta.cl ~ 0.1
          logn.sd = 10
        })
        model({
          ka <- exp(lka + eta.ka)
          cl <- exp(lcl + eta.cl)
          v <- exp(lv + eta.v)
          q <- exp(lq)
          vp <- exp(lvp)
          linCmt() ~ lnorm(logn.sd)
        })
      }

      cmt2fit.logn <- nlmixr(cmt2, dat2, "saem",
                                      control=list(print=0),
                                      table=tableControl(cwres=TRUE, npde=TRUE))

      expect_error(augPred(cmt2fit.logn), NA)

  })

  test_that("augPred with pop only data", {

    one.cmt <- function() {
      ini({
        ## You may label each parameter with a comment
        tka <- 0.45 # Log Ka
        tcl <- log(c(0, 2.7, 100)) # Log Cl
        ## This works with interactive models
        ## You may also label the preceding line with label("label text")
        tv <- 3.45; label("log V")
        ## the label("Label name") works with all models
        add.sd <- 0.7
      })
      model({
        ka <- exp(tka)
        cl <- exp(tcl)
        v <- exp(tv)
        linCmt() ~ add(add.sd)
      })
    }

    fit2 <- nlmixr(one.cmt, nlmixr2data::theo_sd, est="focei",
                   table=tableControl(npde=TRUE))

    expect_error(augPred(fit2), NA)

  })

})

