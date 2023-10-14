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

fit1 <- nlmixr(one.cmt, nlmixr2data::theo_sd, est="nls")


one.cmt <- function() {
  ini({
    tka <- 0.45
    tcl <- log(c(0, 2.7, 100))
    tv <- fix(3.45)
    add.sd <- 0.7
  })
  model({
    ka <- exp(tka)
    cl <- exp(tcl)
    v <- exp(tv)
    linCmt() ~ add(add.sd)
  })
}


one.cmt <- function() {
  ini({
    tka <- 0.45
    tcl <- log(c(0, 2.7, 100))
    tv <- 3.45
    prop.sd <- 0.7
  })
  model({
    ka <- exp(tka)
    cl <- exp(tcl)
    v <- exp(tv)
    linCmt() ~ prop(prop.sd)
  })
}

one.cmt <- function() {
  ini({
    tka <- 0.45
    tcl <- log(c(0, 2.7, 100))
    tv <- 3.45
    prop.sd <- 0.7
    pow.exp <- 1
  })
  model({
    ka <- exp(tka)
    cl <- exp(tcl)
    v <- exp(tv)
    linCmt() ~ pow(prop.sd, pow.exp)
  })
}

one.cmt <- function() {
  ini({
    tka <- 0.45
    tcl <- log(c(0, 2.7, 100))
    tv <- 3.45
    add.sd <- 0.7
    lambda <- c(-2, 0, 2)
  })
  model({
    ka <- exp(tka)
    cl <- exp(tcl)
    v <- exp(tv)
    linCmt() ~ add(add.sd) + yeoJohnson(lambda)
  })
}


#fit1 <- nlmixr(one.cmt, nlmixr2data::theo_sd, est="nlm")
