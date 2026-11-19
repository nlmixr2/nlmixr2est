# The SD-scale IOV fixture of the analytic covariance tests: 24 subjects, two
# doses 48 h apart, and genuine between-occasion variability in CL (sd sqrt(0.08)),
# so iov.cl is identified (theo_sd has none, which drives it onto its boundary).
.covIovData <- function() {
  .testSeed(42)
  .sim <- rxode2::rxode2({
    ka <- exp(0.5 + eta.ka); cl <- exp(1.0 + eta.cl + iov.cl); v <- exp(3.4 + eta.v)
    d/dt(depot)  <- -ka * depot
    d/dt(center) <-  ka * depot - cl / v * center
    cp <- center / v
  })
  .ev <- rxode2::et(seq(0.25, 20, by = 2))
  .ev <- rxode2::et(.ev, seq(48.25, 68, by = 2))
  .ev <- rxode2::et(.ev, amt = 320, time = 0, cmt = "depot")
  .ev <- rxode2::et(.ev, amt = 320, time = 48, cmt = "depot")
  .ev <- as.data.frame(.ev)
  .nid <- 24
  .d <- do.call(
    rbind,
    lapply(seq_len(.nid), function(id) {
      .dd <- .ev
      .dd$iov.cl <- ifelse(.dd$time < 40, rnorm(1, 0, sqrt(0.08)), rnorm(1, 0, sqrt(0.08)))
      .s <- rxode2::rxSolve(
        .sim,
        c(eta.ka = rnorm(1, 0, sqrt(0.3)), eta.cl = rnorm(1, 0, sqrt(0.1)), eta.v = rnorm(1, 0, sqrt(0.1))),
        .dd,
        returnType = "data.frame"
      )
      .obs <- .s[!is.na(.s$cp) & .s$cp > 0, ]
      data.frame(
        ID = id,
        TIME = .obs$time,
        DV = .obs$cp * exp(rnorm(nrow(.obs), 0, 0.1)),
        AMT = 0,
        EVID = 0,
        occ = ifelse(.obs$time < 40, 1L, 2L)
      )
    })
  )
  .dose <- data.frame(
    ID = rep(seq_len(.nid), each = 2),
    TIME = rep(c(0, 48), .nid),
    DV = 0,
    AMT = 320,
    EVID = 1,
    occ = rep(c(1L, 2L), .nid)
  )
  dat <- rbind(.dose, .d)
  dat <- dat[order(dat$ID, dat$TIME, -dat$EVID), ]
  dat
}

# The fixture's model: one IOV eta on CL, additive error
.covIovModel <- function() {
  ini({ tka <- 0.5; tcl <- 1.0; tv <- 3.4
        eta.ka ~ 0.3; eta.cl ~ 0.1; eta.v ~ 0.1
        iov.cl ~ 0.08 | occ
        add.sd <- 0.3 })
  model({ ka <- exp(tka + eta.ka); cl <- exp(tcl + eta.cl + iov.cl); v <- exp(tv + eta.v)
    d/dt(depot)  <- -ka * depot
    d/dt(center) <-  ka * depot - cl / v * center
    cp <- center / v; cp ~ add(add.sd) })
}
