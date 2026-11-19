# Models the BSV correlation tests fit (test-nlmixr2output*.R)
.bsvCorBlockFit <- function() {
  ini({
    tka <- 0.45
    tcl <- 1
    tv <- 3.45
    eta.ka ~ 0.6
    eta.cl + eta.v ~ c(0.3, 0.05, 0.1)
    add.sd <- 0.7
  })
  model({
    ka <- exp(tka + eta.ka)
    cl <- exp(tcl + eta.cl)
    v <- exp(tv + eta.v)
    linCmt() ~ add(add.sd)
  })
}

.bsvCorMixedFit <- function() {
  ini({
    tka <- 0.45
    label("Ka")
    tcl <- 1
    tv <- 3.45
    eta.ka + eta.cl ~ c(0.6, 0.01, 0.3)
    eta.v ~ 0.1
    add.sd <- 0.7
  })
  model({
    ka <- exp(tka) * exp(eta.ka)
    cl <- exp(tcl + eta.cl)
    v <- exp(tv + eta.v)
    linCmt() ~ add(add.sd)
  })
}

.bsvCorIovFit <- function() {
  ini({
    tka <- 0.45
    tcl <- 1
    tv <- 3.45
    add.sd <- 0.7
    eta.cl + eta.v ~ c(0.3, 0.05, 0.1)
    iov.cl + iov.v ~ c(0.1, 0.01, 0.1) | occ
  })
  model({
    ka <- exp(tka)
    cl <- exp(tcl + eta.cl + iov.cl)
    v <- exp(tv + eta.v + iov.v)
    linCmt() ~ add(add.sd)
  })
}

# tcl is mu-referenced by two etas, so neither eta is named by it
.bsvCorSharedFit <- function() {
  ini({
    tcl <- 1
    tv <- 3.45
    eta.cl + eta.cl2 + eta.v ~ c(0.3, 0.02, 0.2, 0.03, 0.04, 0.1)
    add.sd <- 0.7
  })
  model({
    cl <- exp(tcl + eta.cl)
    ka <- exp(tcl + eta.cl2)
    v <- exp(tv + eta.v)
    linCmt() ~ add(add.sd)
  })
}

.bsvCorSameFit <- function() {
  ini({
    tka <- 0.45
    tcl <- 1
    tv <- 3.45
    tq <- 1
    tvp <- 3
    add.sd <- 0.7
    eta.cl + eta.v ~ c(0.3, 0.05, 0.1)
    eta.q + eta.vp ~ same()
  })
  model({
    ka <- exp(tka)
    cl <- exp(tcl + eta.cl)
    v <- exp(tv + eta.v)
    q <- exp(tq + eta.q)
    vp <- exp(tvp + eta.vp)
    linCmt() ~ add(add.sd)
  })
}

# saem refuses correlated occasion etas
.bsvCorIovSaemFit <- function() {
  ini({
    tka <- 0.45
    tcl <- 1
    tv <- 3.45
    add.sd <- 0.7
    eta.cl + eta.v ~ c(0.3, 0.05, 0.1)
    iov.cl ~ 0.1 | occ
  })
  model({
    ka <- exp(tka)
    cl <- exp(tcl + eta.cl + iov.cl)
    v <- exp(tv + eta.v)
    linCmt() ~ add(add.sd)
  })
}
