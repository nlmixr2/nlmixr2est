model.2cpt.ode.wtcl.sexv2 <- function() {
  ini({
    tka <- log(1.14)
    tcl <- log(0.0190)
    tv2  <- log(2.12)
    tv3  <- log(20.4)
    tq   <- log(0.383)
    wteff  <- 0.35
    sexeff <- -0.2
    eta.ka ~ 0.1
    eta.cl ~ 0.2
    eta.v2 ~ 0.3
    eta.v3 ~ 0.4
    eta.q ~ 0.5
    prop.err <- 0.075
  })
  model({
    ka <- exp(tka + eta.ka)
    cl <- exp(tcl + wteff*log(WT/70) + eta.cl)
    v2 <- exp(tv2 + sexeff*(SEX) + eta.v2)
    v3 <- exp(tv3 + eta.v3)
    q  <- exp(tq + eta.q)
    d/dt(depot) = -ka * depot
    d/dt(center) = ka * depot - cl / v2 * center + q/v3 * periph - q/v2 * center
    d/dt(periph) = q/v2 * center - q/v3 * periph
    cp = center / v2
    cp ~ prop(prop.err)
  })
}

f <- nlmixr(model.2cpt.ode.wtcl.sexv2)

 model.2cpt.ode.wtcl.sexv2 <- function() {
  ini({
    tka <- log(1.14)
    tcl <- log(0.0190)
    tv2  <- log(2.12)
    tv3  <- log(20.4)
    tq   <- log(0.383)
    wteff  <- 0.35
    sexeff <- -0.2
    eta.ka + eta.cl ~ c(0.1,
                        0.001, 0.2)
    eta.v2 ~ 0.3
    eta.v3 ~ 0.4
    eta.q ~ 0.5
    prop.err <- 0.075
  })
  model({
    ka <- exp(tka + eta.ka)
    cl <- exp(tcl + wteff*log(WT/70) + eta.cl)
    v2 <- exp(tv2 + sexeff*(SEX) + eta.v2)
    v3 <- exp(tv3 + eta.v3)
    q  <- exp(tq + eta.q)
    d/dt(depot) = -ka * depot
    d/dt(center) = ka * depot - cl / v2 * center + q/v3 * periph - q/v2 * center
    d/dt(periph) = q/v2 * center - q/v3 * periph
    cp = center / v2
    cp ~ prop(prop.err)
  })
}

f <- nlmixr(model.2cpt.ode.wtcl.sexv2)




model.2cpt.cf.wtcl.sexv2 <- function() {
  ini({
    tka <- log(1.5)
    tcl <- log(0.05)
    tv2  <- log(2.5)
    tv3  <- log(25)
    tq   <- log(0.5)
    wteff <- 0.01
    sexeff <- -0.01
    eta.ka ~ 1
    eta.cl ~ 1
    eta.v2 ~ 1
    eta.v3 ~ 1
    eta.q ~ 1
    prop.err <- 0.5
  })
  model({
    ka <- exp(tka + eta.ka)
    cl <- exp(tcl + wteff*lwt70 + eta.cl)
    v2 <- exp(tv2 + sexeff*SEX + eta.v2)
    v3 <- exp(tv3 + eta.v3)
    q  <- exp(tq + eta.q)
    linCmt() ~ prop(prop.err)
  })
}

f <- nlmixr(model.2cpt.cf.wtcl.sexv2)
## > f$nlme.specs
## $fixed
## tka + tcl + tv2 + tv3 + tq + wteff + sexeff ~ 1
## <environment: 0x5642b95d2638>

## $random
## Positive definite matrix structure of class pdDiag representing
##      [,1] [,2] [,3] [,4] [,5]
## [1,]    1    0    0    0    0
## [2,]    0    1    0    0    0
## [3,]    0    0    1    0    0
## [4,]    0    0    0    1    0
## [5,]    0    0    0    0    1

## $start
##        tka        tcl        tv2        tv3         tq      wteff     sexeff
##  0.4054651 -2.9957323  0.9162907  3.2188758 -0.6931472  0.0100000 -0.0100000

## > `
