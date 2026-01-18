nmTest({
  ## define model

  my.model.pk <- function() {
    ini({
      lV <- 0.795860222982066
      lka <- -2.48182829010656
      ln <- 3.10678040444175
      lCL <- 1.70924362519995
      lktr <- 4.37261675943197
      err.pk <- c(0, 0.949673109350856)
      eta.lka ~ 0.0195139946086797
      eta.lCL ~ 0.0446009979019879
      eta.ln ~ 0.81033197903365
    })
    model({
      V <- exp(lV)
      CL <- exp(lCL + eta.lCL)
      ka <- exp(lka + eta.lka)
      ktr <- exp(lktr)
      n <- exp(ln + eta.ln)
      bio <- 1
      d/dt(depot) <- exp(log(bio * podo(depot)) + log(ktr) +
                           n * log(ktr * tad(depot)) - ktr * tad(depot) - lgamma(n +
                                                                                   1)) - ka * depot
      d/dt(central) <- ka * depot - CL/V * central
      central_conc <- central/V
      central_conc ~ add(err.pk)
    })
  }
  ## fit the model
  d <- theo_md |>
    dplyr::mutate(EVID=ifelse(EVID == 0, 0L, 7L))

  tmp <- expect_error(nlmixr(my.model.pk, d, est="focei", control=foceiControl(print=0L)), NA)

  ## Now try with saem
  tmp <- expect_error(nlmixr(my.model.pk, d, est="saem", control=saemControl(print=0L)), NA)
})
