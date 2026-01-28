nmTest({
  test_that("test augPred", {

    PKdata <- nlmixr2data::warfarin |>
      dplyr::filter(dvid == "cp") |>
      dplyr::select(-dvid) |>
      dplyr::mutate(sex = ifelse(sex == "male", 1, 0))

    One.comp.KA.solved <- function() {
      ini({
        lka  <- log(1.15)
        lcl  <- log(0.135)
        lv   <- log(8)
        prop.err <- 0.15
        add.err  <- 0.6
        eta.ka ~ 0.5
        eta.cl ~ 0.1
        eta.v  ~ 0.1
      })
      model({
        cl <- exp(lcl + eta.cl)
        v  <- exp(lv + eta.v)
        ka <- exp(lka + eta.ka)
        linCmt() ~ prop(prop.err) + add(add.err)
      })
    }

    fitOne.comp.KA.solved_S <-
      .nlmixr(
        One.comp.KA.solved,
        PKdata,
        est = "saem",
        control = saemControlFast,
        tableControl(cwres = TRUE, npde=TRUE)
      )
    expect_error(augPred(fitOne.comp.KA.solved_S), NA)

    ap <- augPred(fitOne.comp.KA.solved_S)

    expect_equal(as.character(ap[ap$id == 1 & ap$time == 120, "ind"]),
                 c("Individual", "Population"))

      df <-
        tibble::tibble(
          ID = c(rep(1, 6), rep(2, 6)),
          TIME = c(0.00, 12.11, 18.41, 23.89, 36.00, 43.51, 0.00, 12.00, 20.00, 24.00, 36.80, 45.00),
          AMT = c(1000, 1000, NA, 1000, 1000, NA, 1000, 2000, NA, 1000, 1000, NA),
          DUR = c(2.5, 2.5, NA, 2.5, 2.5, NA, 2.5, 2.5, NA, 2.5, 2.5, NA),
          DV = c(NA, NA, 3.0, NA, NA, 9.6, NA, NA, 7.0, NA, NA, 2.8),
          WT = c(rep(55, 6), rep(48, 6))
        ) |>
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

      fit <- .nlmixr(fun, df, list(print=0), est="posthoc")

      expect_error(augPred(fit), NA)
  })

  test_that("test augPred with xgxr dataset", {

    dat <- xgxr::case1_pkpd |>
      dplyr::rename(DV=LIDV) |>
      dplyr::filter(CMT %fin% 1:2) |>
      dplyr::filter(TRTACT != "Placebo")

      doses <- unique(dat$DOSE)
      nid <- 3 # 7 ids per dose group
      dat2 <- do.call("rbind",
                      lapply(doses, function(x) {
                        ids <- dat |>
                          dplyr::filter(DOSE == x) |>
                          dplyr::reframe(ids=unique(ID)) |>
                          dplyr::pull()
                        ids <- ids[seq(1, nid)]
                        dat |>
                          dplyr::filter(ID %fin% ids)
                      }))

      # Use centralized model from helper-models.R
      cmt2 <- two.compartment

      cmt2fit.logn <-
        .nlmixr(
          cmt2, dat2, "saem",
          control = saemControlFast,
          table=tableControl(cwres=TRUE, npde=TRUE)
        )

      expect_error(augPred(cmt2fit.logn), NA)
  })

  test_that("augPred with pop only data", {
    one.cmt <- function() {
      ini({
        tka <- 0.45
        tcl <- log(c(0, 2.7, 100))
        tv <- 3.45
        add.sd <- 0.7
      })
      model({
        ka <- exp(tka)
        cl <- exp(tcl)
        v <- exp(tv)
        linCmt() ~ add(add.sd)
      })
    }

    fit2 <-
      .nlmixr(
        one.cmt, nlmixr2data::theo_sd, est="focei",
        control = foceiControl(eval.max = 1),
        table=tableControl(npde=TRUE)
      )
    expect_error(augPred(fit2), NA)
  })

  test_that("mixed pkpd with effect compartment augpred", {

    dat <- nlmixr2data::warfarin

    mod <- function () {
      ini({
        tktr <- -0.0407039444259225
        tcl <- -1.94598426244892
        tv <- 2.09185800199064
        eps.pkprop <- c(0, 0.103324277742915)
        eps.pkadd <- c(0, 0.41909330805883)
        tc50 <- 0.121920646717701
        tkout <- -3.1790123282253
        te0 <- 4.43865746919035
        eps.pdadd <- c(0, 5.96046447753906e-07)
        eta.ktr ~ 0.625781701127507
        eta.cl ~ 0.102936117458982
        eta.v ~ 0.0491541297805943
        eta.c50 ~ 0.41589473728507
        eta.kout ~ 0.108340169259417
        eta.e0 ~ 0.0632084214464694
      })
      model({
        ktr <- exp(tktr + eta.ktr)
        cl <- exp(tcl + eta.cl)
        v <- exp(tv + eta.v)
        c50 <- exp(tc50 + eta.c50)
        kout <- exp(tkout + eta.kout)
        e0 <- exp(te0 + eta.e0)
        emax <- 1
        cp <- central/v
        d/dt(depot) <- -ktr * depot
        d/dt(central) <- ktr * trans - cl * cp
        d/dt(trans) <- ktr * depot - ktr * trans
        d/dt(ce) = kout * (cp - ce)
        effect <- e0 * (1 - emax * ce/(c50 + ce))
        cp ~ prop(eps.pkprop) + add(eps.pkadd) | cp
        effect ~ add(eps.pdadd) | pca
      })
    }

    fit <- .nlmixr(mod, dat, "posthoc")

    expect_error(augPred(fit), NA)
  })

  test_that("augPred with zero etas", {

    # Use centralized model from helper-models.R
    model.1compt.depot1 <- one.compartment |>
      ini(eta.ka~0)

    fit1 <- .nlmixr(model.1compt.depot1, theo_sd, est="focei",
                    foceiControlFast)

    expect_error(augPred(fit1), NA)

  })
})
