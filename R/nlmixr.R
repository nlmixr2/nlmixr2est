
##' @rdname nlmixr2Est000
##'@export
nlmixr2Est000.nlme <- function(env, ...) {
  with(env, {
    if (.nid <= 1) stop("nlme is for mixed effects models, try 'dynmodel' (need more than 1 individual)")
    if (.nTv != 0) stop("nlme does not support time-varying covariates (yet)")
    data <- .as.data.frame(data)
    if (length(uif$predDf$cond) > 1) stop("nlmixr2 nlme does not support multiple endpoints.")
    pt <- proc.time()
    est.type <- est
    if (est == "nlme.free") {
      fun <- uif$nlme.fun
      specs <- uif$nlme.specs
    } else if (est == "nlme.mu") {
      fun <- uif$nlme.fun.mu
      specs <- uif$nlme.specs.mu
    } else if (est == "nlme.mu.cov") {
      fun <- uif$nlme.fun.mu.cov
      specs <- uif$nlme.specs.mu.cov
    } else {
      if (!is.null(uif$nlme.fun.mu.cov)) {
        est.type <- "nlme.mu.cov"
        fun <- uif$nlme.fun.mu.cov
        specs <- uif$nlme.specs.mu.cov
      } else if (!is.null(uif$nlme.fun.mu)) {
        est.type <- "nlme.mu"
        fun <- uif$nlme.fun.mu
        specs <- uif$nlme.specs.mu
      } else {
        est.type <- "nlme.free"
        fun <- uif$nlme.fun
        specs <- uif$nlme.fun.specs
      }
    }
    grp.fn <- uif$grp.fn
    dat$nlmixr2.grp <-
      factor(apply(dat, 1, function(x) {
        cur <- x
        names(cur) <- names(dat)
        with(as.list(cur), {
          return(grp.fn())
        })
      }))
    dat$nlmixr2.num <- seq_along(dat$nlmixr2.grp)
    .addProp <- "combined2"
    if (!is.null(control$addProp)) .addProp <- control$addProp
    if (!any(.addProp == c("combined2", "combined1"))) stop("addProp needs to either be 'combined1' and 'combined2'")
    uif$env$.addProp <- .addProp
    weight <- uif$nlme.var
    if (sum.prod) {
      rxode <- rxode2::rxSumProdModel(uif$rxode.pred)
    } else {
      rxode <- uif$rxode.pred
    }
    .atol <- 1e-8
    if (!is.null(control$atol)) .atol <- control$atol
    .rtol <- 1e-8
    if (!is.null(control$rtol)) .rtol <- control$rtol
    .maxsteps <- 5000
    if (!is.null(control$maxstepsOde)) .maxsteps <- control$maxstepsOde
    if (is(weight, "varConstProp")) {
      control$sigma <- 1
    }
    fit <- nlme_ode(dat,
                    model = rxode,
                    par_model = specs,
                    par_trans = fun,
                    response = "nlmixr2_pred",
                    weight = weight,
                    verbose = TRUE,
                    control = control,
                    atol = .atol,
                    rtol = .rtol,
                    maxsteps = .maxsteps,
                    ...
                    )
    class(fit) <- c(est.type, class(fit))
    .ret <- try({
      as.focei.nlmixr2Nlme(fit, uif, pt, data = dat, calcResid = calc.resid, nobs2 = nobs2,
                          keep=.keep, drop=.drop, IDlabel=.lab, table=table)
    })
    if (inherits(.ret, "try-error")) {
      warning("Error converting to nlmixr2 UI object, returning nlme object")
      return(fit)
    }
    if (inherits(.ret, "nlmixr2FitCore")) {
      .ret <- nlmixr2FitUpdateParams(.ret, origData = .origData)
    }
    if (inherits(.ret, "nlmixr2FitCore")) {
      .env <- .ret$env
      assign("startTime", start.time, .env)
      assign("est", est, .env)
      assign("stopTime", Sys.time(), .env)
      assign("origControl", control, .env)
      assign("modelId", .modelId, .env)
    }
    return(.ret)
  })
}
