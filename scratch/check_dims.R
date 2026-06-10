library(devtools)
load_all()

set.seed(42)
n_subj <- 30
sub_pop <- rbinom(n_subj, 1, 0.6) + 1 # 1 or 2 (40% sub-pop 1, 60% sub-pop 2)
cl_sim <- ifelse(sub_pop == 1, 1.2, 6.0)

sim_data <- do.call(rbind, lapply(1:n_subj, function(i) {
  subj_cl <- cl_sim[i]
  times <- c(0.5, 1, 2, 4, 8, 12, 24)
  ka_val <- 1.5
  v_val <- 24.0
  k_val <- subj_cl / v_val
  cp <- 100 * ka_val / (v_val * (ka_val - k_val)) * (exp(-k_val * times) - exp(-ka_val * times)) + rnorm(length(times), 0, 0.05)
  cp[cp < 0] <- 0
  data.frame(
    ID = i,
    TIME = c(0, times),
    AMT = c(100, rep(0, length(times))),
    EVID = c(1, rep(0, length(times))),
    DV = c(0, cp),
    CMT = c(1, rep(2, length(times)))
  )
}))

one.compartment.mix <- function() {
  ini({
    tka <- log(1.5)
    tcl1 <- log(1.0)
    tcl2 <- log(5.0)
    tv <- log(20)
    p1 <- 0.5
    eta.cl ~ 0.01
    eta.v ~ 0.01
    eta.ka ~ 0.01
    add.sd <- 0.05
  })
  model({
    ka <- exp(tka + eta.ka)
    cl <- mix(exp(tcl1 + eta.cl), p1, exp(tcl2 + eta.cl))
    v <- exp(tv + eta.v)
    d/dt(depot) <- -ka * depot
    d/dt(center) <- ka * depot - cl / v * center
    cp <- center / v
    cp ~ add(add.sd)
  })
}

ui <- rxode2::as.rxUi(one.compartment.mix)
ui <- nlmixr2est:::.rxUiDecompressModelFun(ui)

.model <- ui$saemModelList
.inits <- ui$saemInit
print(ui$saemModelOmegaFixed)
print(ui$saemModelOmegaFixedValues)
.rxControl <- rxode2::rxControl()
.ue <- nlmixr2est:::.uninformativeEtas(ui, handleUninformativeEtas=TRUE, data=sim_data, attr(.model$saem_mod, "rx"), rxControl=.rxControl)

cfg <- nlmixr2est:::.configsaem(model=.model,
                               data=sim_data,
                               inits=.inits,
                               mcmc=list(niter = c(200, 300), nmc = 3, nu = c(2, 2, 2)),
                               rxControl=.rxControl,
                               distribution="normal",
                               fixedOmega=ui$saemModelOmegaFixed,
                               fixedOmegaValues=ui$saemModelOmegaFixedValues,
                               parHistThetaKeep=ui$saemParHistThetaKeep,
                               parHistOmegaKeep=ui$saemParHistOmegaKeep,
                               seed=99,
                               DEBUG=0,
                               tol=1e-6,
                               itmax=30,
                               type="nelder-mead",
                               lambdaRange=3,
                               powRange=10,
                               odeRecalcFactor=10^0.5,
                               maxOdeRecalc=10^0.5,
                               indTolRelax=TRUE,
                               nres=ui$saemModNumEst,
                               perSa=0.75,
                               perNoCor=0.75,
                               perFixOmega=0.1,
                               perFixResid=0.1,
                               resFixed=ui$saemResFixed,
                               ue=.ue,
                               mixProb=rxUiGet.saemMixProb(list(ui)))

cat("cfg$nlambda1:", cfg$nlambda1, "\n")
cat("cfg$nlambda0:", cfg$nlambda0, "\n")
cat("cfg$ilambda1:", paste(cfg$ilambda1, collapse=", "), "\n")
cat("cfg$ilambda0:", paste(cfg$ilambda0, collapse=", "), "\n")
cat("ui$saemParHistThetaKeep:", paste(ui$saemParHistThetaKeep, collapse=", "), "\n")
cat("sum(ui$saemParHistThetaKeep):", sum(ui$saemParHistThetaKeep), "\n")
cat("ncol(cfg$par.hist):", ncol(cfg$par.hist), "\n")

