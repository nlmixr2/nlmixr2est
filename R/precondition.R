.getUiFunFromIniAndModel <- function(ui, ini, model) {
  .ls <- ls(ui$meta, all.names=TRUE)
  .ret <- vector("list", length(.ls) + 3)
  .ret[[1]] <- quote(`{`)
  for (.i in seq_along(.ls)) {
    .ret[[.i + 1]] <- eval(parse(text=paste("quote(", .ls[.i], "<-", deparse1(ui$meta[[.ls[.i]]]), ")")))
  }
  .len <- length(.ls)
  .ret[[.len + 2]] <- ini
  .ret[[.len + 3]] <- model
  .retf <- function(){}
  body(.retf) <- as.call(.ret)
  .retf
}

#' Linearly re-parameterize the model to be less sensitive to rounding errors
#'
#' @param fit A nlmixr2 fit to be preconditioned
#' @param estType Once the fit has been linearly reparameterized,
#'   should a "full" estimation, "posthoc" estimation or simply a
#'   estimation of the covariance matrix "none" before the fit is
#'   updated
#' @param ntry number of tries before giving up on a pre-conditioned
#'   covariance estimate
#'
#' @return A nlmixr2 fit object that was preconditioned to stabilize
#'   the variance/covariance calculation
#'
#' @export
#'
#' @references Aoki Y, Nordgren R, Hooker AC. Preconditioning of
#'   Nonlinear Mixed Effects Models for Stabilisation of
#'   Variance-Covariance Matrix Computations. AAPS
#'   J. 2016;18(2):505-518. doi:10.1208/s12248-016-9866-5
#'
preconditionFit <- function(fit, estType = c("full", "posthoc", "none"),
                            ntry = 10L) {
  rxode2::.setWarnIdSort(FALSE)
  on.exit(rxode2::.setWarnIdSort(TRUE))
  if (!exists("R", fit$env)) {
    stop("this assumes a covariance matrix with a R matrix",
      call. = FALSE
    )
  }
  .R <- fit$R
  .covMethod <- ""
  .i <- 1
  while (.i < ntry & .covMethod != "r,s") {
    .i <- .i + 1
    pre <- preCondInv(.R)
    P <- symengine::Matrix(pre)
    d0 <- dimnames(fit$R)[[1]]
    d <- paste0("nlmixr2Pre_", dimnames(fit$R)[[1]])
    d2 <- symengine::Matrix(d)
    modExtra <- paste(paste0(d0, "=", sapply(as.vector(P %*% d2), as.character)), collapse = "\n")
    preInv <- solve(pre)

    .ini <- as.data.frame(fit$ui$iniDf)
    newEst <- setNames(as.vector(preInv %*% matrix(fit$theta[d0])), d0)
    for (v in d0) {
      .w <- which(.ini$name == v)
      .ini$lower[.w] <- -Inf
      .ini$upper[.w] <- Inf
      .ini$est[.w] <- newEst[v]
      .ini$name[.w] <- paste0("nlmixr2Pre_", v)
    }
    .ini <- as.expression(lotri::as.lotri(.ini))
    .ini[[1]] <- quote(`ini`)
    .newModel <- eval(parse(text = paste0("quote(model({", modExtra, "\n", fit$ui$fun.txt, "}))")))
    .newModel <- .getUiFunFromIniAndModel(fit$ui, .ini, .newModel)
    .newModel <- .newModel()
    .ctl <- fit$control
    estType <- match.arg(estType)
    if (estType == "none") {
      .ctl$maxInnerIterations <- 0
      .ctl$maxOuterIterations <- 0
      .ctl$boundTol <- 0
      .ctl$etaMat <- as.matrix(fit$eta[, -1])
      .ctl$calcTables <- FALSE
      .ctl$compress <- FALSE
    } else if (estType == "posthoc") {
      .ctl$maxOuterIterations <- 0
      .ctl$boundTol <- 0
      .ctl$calcTables <- FALSE
      .ctl$compress <- FALSE
    } else if (estType == "full") {
      .ctl$boundTol <- 0
      .ctl$calcTables <- FALSE
      .ctl$compress <- FALSE
    }
    .ctl$covMethod <- 1L
    ## FIXME compare objective functions
    newFit <- suppressWarnings(nlmixr2(.newModel, getData(fit), est = "focei", control = .ctl))
    .R <- newFit$R
    .covMethod <- newFit$covMethod
  }
  if (.covMethod != "r,s") {
    stop("preconditioning failed after ", ntry, "tries",
      call. = FALSE
    )
  }
  cov <- pre %*% newFit$cov %*% t(pre)
  dimnames(cov) <- dimnames(pre)
  assign("precondition", cov, envir = fit$env)
  .setCov(fit, covMethod = cov)
  assign("covMethod", "precondition", envir=fit$env)
  return(invisible(fit$env$precondition))
}
