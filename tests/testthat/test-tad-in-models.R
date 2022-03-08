skip_on_cran()
## define model
my.model.pk <- function(){
  ini({
    lV   <- 0.2
    lka  <- -2
    ln   <- 0.3
    lCL  <- 0.9
    lktr <- 3
    eta.lka  ~ 0.5
    eta.lCL  ~ 0.5
    eta.ln   ~ 0.5
    err.pk <- 0.1
  })

  model({
    V   <- exp(lV)
    CL  <- exp(lCL + eta.lCL)
    ka  <- exp(lka + eta.lka)
    ktr <- exp(lktr)
    n   <- exp(ln + eta.ln)
    # mtt <- (n+1)/ktr
    bio <- 1
    d/dt(depot)   <-  exp(log(bio*podo(depot))+log(ktr)+n*log(ktr*tad(depot))-ktr*tad(depot)-lgamma(n+1)) - ka * depot
    d/dt(central) <-  ka * depot - CL/V * central
    central_conc <- central / V
    central_conc ~ add(err.pk)
  })
}

## fit the model
d <- theo_md %>%
  dplyr::mutate(EVID=ifelse(EVID == 0, 0L, 7L))

tmp <- expect_error(nlmixr(my.model.pk, theo_md, est="focei", control=foceiControl(print=0L)), NA)

## Now try with saem
tmp <- expect_error(nlmixr(my.model.pk, theo_md, est="saem", control=saemControl(print=0L)), NA)
