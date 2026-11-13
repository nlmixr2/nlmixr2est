# Wall time of a SAEM fit at a given thread count, with a fingerprint of the
# estimates so builds can be checked for bit-identical results (developer
# script; not run on CI/CRAN).  See README.md in this folder.
#
# Usage: Rscript design/saem-readback/bench.R <label> <cores> [nrep] [nsub] [nBurn] [nEm]
suppressMessages(library(nlmixr2est))
args <- commandArgs(TRUE)
label <- args[1]
cores <- as.integer(args[2])
nrep <- if (length(args) > 2) as.integer(args[3]) else 3L
nsub <- if (length(args) > 3) as.integer(args[4]) else 200L
nBurn <- if (length(args) > 4) as.integer(args[5]) else 60L
nEm <- if (length(args) > 5) as.integer(args[6]) else 60L

mod <- function() {
  ini({
    tka <- log(1.1)
    tcl <- log(3)
    tv1 <- log(30)
    tq <- log(5)
    tv2 <- log(60)
    eta.ka ~ 0.1
    eta.cl ~ 0.09
    eta.v1 ~ 0.06
    prop.sd <- 0.1
  })
  model({
    ka <- exp(tka + eta.ka)
    cl <- exp(tcl + eta.cl)
    v1 <- exp(tv1 + eta.v1)
    q <- exp(tq)
    v2 <- exp(tv2)
    d / dt(depot) <- -ka * depot
    d / dt(central) <- ka * depot - (cl + q) / v1 * central + q / v2 * periph
    d / dt(periph) <- q / v1 * central - q / v2 * periph
    cp <- central / v1
    cp ~ prop(prop.sd)
  })
}

ev <- rxode2::et(amt = 100, ii = 24, addl = 2)
ev <- rxode2::et(ev, c(0.25, 0.5, 1, 2, 4, 8, 12, 24, 48, 49, 52, 72))
sim <- rxode2::rxWithSeed(
  2026,
  rxode2::rxSolve(
    mod,
    ev,
    nSub = nsub,
    returnType = "data.frame",
    addDosing = TRUE,
    cores = 1L
  )
)
idCol <- if ("sim.id" %in% names(sim)) "sim.id" else "id"
dat <- data.frame(
  ID = sim[[idCol]],
  TIME = sim$time,
  EVID = ifelse(is.na(sim$evid), 0, sim$evid),
  AMT = ifelse(is.na(sim$amt), 0, sim$amt),
  DV = ifelse(is.na(sim$evid) | sim$evid == 0, sim$sim, NA)
)
dat$EVID[dat$AMT > 0] <- 1

fingerprint <- function(x) {
  f <- tempfile()
  on.exit(unlink(f))
  writeLines(paste(sprintf("%a", x), collapse = ","), f)
  substr(unname(tools::md5sum(f)), 1, 12)
}

times <- numeric(nrep)
prints <- character(nrep)
for (r in seq_len(nrep)) {
  t0 <- proc.time()[["elapsed"]]
  fit <- suppressMessages(suppressWarnings(
    nlmixr2(
      mod,
      dat,
      est = "saem",
      control = saemControl(
        print = 0,
        nBurn = nBurn,
        nEm = nEm,
        calcTables = FALSE,
        covMethod = "",
        logLik = FALSE,
        seed = 99,
        rxControl = rxode2::rxControl(cores = cores)
      )
    )
  ))
  times[r] <- proc.time()[["elapsed"]] - t0
  prints[r] <- fingerprint(c(fit$theta, as.vector(fit$omega), unlist(fit$eta[, -1])))
}
cat(sprintf(
  "%s cores=%d: median %.3f s (min %.3f max %.3f) fp %s\n",
  label,
  cores,
  median(times),
  min(times),
  max(times),
  paste(unique(prints), collapse = " ")
))
