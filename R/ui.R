#' @importFrom rxode2 ini
#' @export
ini <- rxode2::ini

#' @importFrom rxode2 model
#' @export
model <- rxode2::model

#' @importFrom rxode2 lotri
#' @export
lotri <- rxode2::lotri


nlmixr2findRhs <- function(x) {
  ## Modified from http://adv-r.had.co.nz/Expressions.html find_assign4
  if (is.atomic(x)) {
    character()
  } else if (is.name(x)) {
    return(as.character(x))
  } else if (is.call(x)) {
    if ((identical(x[[1]], quote(`<-`)) ||
      identical(x[[1]], quote(`=`))) &&
        is.name(x[[2]])) {
      unique(c(unlist(nlmixr2findRhs(x[[3]]))))
    } else {
      x1 <- x[-1]
      unique(unlist(lapply(x1, nlmixr2findRhs)))
    }
  } else if (is.pairlist(x)) {
    unique(unlist(lapply(x, nlmixr2findRhs)))
  } else {
    stop("Don't know how to handle type ", typeof(x),
         call. = FALSE
         )
  }
}

nlmixr2findLhs <- function(x) {
  ## Modified from http://adv-r.had.co.nz/Expressions.html find_assign4
  if (is.atomic(x) || is.name(x)) {
    character()
  } else if (is.call(x)) {
    if ((identical(x[[1]], quote(`<-`)) ||
      identical(x[[1]], quote(`=`))) &&
      is.name(x[[2]])) {
      lhs <- as.character(x[[2]])
    } else {
      lhs <- character()
    }
    unique(c(lhs, unlist(lapply(x, nlmixr2findLhs))))
  } else if (is.pairlist(x)) {
    unique(unlist(lapply(x, nlmixr2findLhs)))
  } else {
    stop("Don't know how to handle type ", typeof(x),
      call. = FALSE
    )
  }
}

nlmixr2findRhsLhs <- function(x) {
  .lhs <- nlmixr2findLhs(x)
  .rhs <- nlmixr2findRhs(x)
  .rhs <- setdiff(.rhs, .lhs)
  list(lhs=.lhs, rhs=.rhs)
}

.deparse <- function(expr) {
  deparse(expr, width.cutoff = 500, control = "useSource")
}

.bodyDewrap <- function(ret) {
  .ret <- ret
  if (length(.ret) > 1) {
    .ret[1] <- sub(rex::rex(
      start, or(group("function", any_spaces, "(", any_spaces, ")", any_spaces), ""),
      any_spaces, any_of("("), any_spaces, "{", any_spaces
    ), "", .ret[1], perl = TRUE)
    if (.ret[1] == "") .ret <- .ret[-1]
  }
  if (length(.ret) > 1) {
    .len <- length(.ret)
    .ret[.len] <- sub(rex::rex(any_spaces, "}", any_spaces, any_of(")"), any_spaces, end), "", .ret[.len])
    if (.ret[.len] == "") .ret <- .ret[-.len]
  }
  return(.ret)
}

.deparse1 <- function(expr) {
  return(.bodyDewrap(deparse(expr, width.cutoff = 500)))
}

.pipedData <- NULL
.clearPipedData <- function() {
  assignInMyNamespace(".pipedData", NULL)
}

.getPipedData <- function() {
  .pipedData
}

.getUif <- function(what) {
  if (inherits(what, "nlmixr2FitCore")) {
    .uif <- what$uif
    assignInMyNamespace(".pipedData", getData(what))
  } else if (inherits(what, "nlmixr2UI")) {
    .uif <- what
  } else if (inherits(what, "function")) {
    .uif <- nlmixr2UI(what)
  } else {
    stop("Do not know how to handle object")
  }
  return(.uif)
}

.processLotri <- function(.code, .uif) {
  .mat <- try(eval(parse(text = sprintf("lotri::lotri(%s)", .code))), silent = TRUE)
  if (inherits(.mat, "matrix")) {
    .d <- dimnames(.mat)[[1]]
    .ini2 <- .as.data.frame(.uif$ini)
    .ini1 <- .ini2[is.na(.ini2$neta1), ]
    .ini2 <- .ini2[!is.na(.ini2$neta1), ]
    .ini4 <- .ini2[!is.na(.ini2$err), ]
    .ini2 <- .ini2[is.na(.ini2$err), ]
    .m2 <- as.vector(.mat)
    .l2 <- length(.d) * 2
    .df <- .data.frame(n1 = rep(.d, each = length(.d)), n2 = rep(.d, length(.d)), val = .m2)
    .dfI <- .as.data.frame(.uif$ini[, c("neta1", "neta2", "name")])
    .dfI <- .dfI[!is.na(.dfI$neta1), ]
    .dfI <- .dfI[which(.dfI$neta1 == .dfI$neta2), c("neta1", "name")]
    .d2 <- paste(.dfI$name)
    .diff <- setdiff(.d, .d2)
    if (length(.diff) > 0) {
      stop(sprintf(
        "trying to provide an estimate for non-existant eta: %s",
        paste(.diff, collapse = ", ")
      ))
    }
    names(.dfI)[2] <- "n1"
    .df$n1 <- paste(.df$n1)
    .dfI$n1 <- paste(.dfI$n1)
    .df <- merge(.df, .dfI, all.x = TRUE)
    names(.dfI) <- c("neta2", "n2")
    .df <- merge(.df, .dfI, all.x = TRUE)
    .df <- .df[order(.df$neta1, .df$neta2), ]
    ## Name becomes (eta.cl,eta.ka)
    .df$name <- ifelse(.df$neta1 == .df$neta2, .df$n1, paste0("(", .df$n1, ",", .df$n2, ")"))
    .df$name2 <- ifelse(.df$neta1 == .df$neta2, .df$n1, paste0("(", .df$n2, ",", .df$n1, ")"))
    .ini3 <- do.call("rbind", lapply(seq_along(.df$name), function(i) {
      .name1 <- paste(.df$name[i])
      .name2 <- paste(.df$name2[i])
      .w <- which(.name1 == .ini2$name)
      if (length(.w) == 1) {
        .ini2$est[.w] <<- .df$val[i]
        return(NULL)
      }
      .w <- which(.ini2$name == .name2)
      if (length(.w) == 1) {
        .ini2$est[.w] <<- .df$val[i]
        return(NULL)
      }
      if (.df$neta1[i] < .df$neta2[i]) {
        return(NULL)
      }
      .df2 <- nlmixr2BoundsTemplate
      .df2$neta1 <- .df$neta1[i]
      .df2$neta2 <- .df$neta2[i]
      .df2$name <- .df$name[i]
      .df2$est <- .df$val[i]
      .df2$lower <- -Inf
      .df2$upper <- Inf
      .df2$fix <- FALSE
      .df2$condition <- "ID"
      return(.df2)
    }))
    .ini2 <- rbind(.ini2, .ini3)
    .ini2 <- .ini2[order(.ini2$neta1, .ini2$neta2), ]
    .ini2 <- rbind(.ini1, .ini2, .ini4)
    class(.ini2) <- c("nlmixr2Bounds", "data.frame")
    .uif$ini <- .ini2
    return(.uif)
  } else {
    stop(sprintf("invalid syntax: %s", .code))
  }
}

#'
#'
#'
##' Create the nlme specs list for nlmixr2 nlme solving
##' @inheritParams nlmixr2UI.nlmefun
##' @param mu.type is the mu-referencing type of model hat nlme will be using.
##' @return specs list for nlme
##' @author Matthew L. Fidler
nlmixr2UI.nlme.specs <- function(object, mu.type = c("thetas", "covariates", "none")) {
  mu.type <- match.arg(mu.type)
  if (mu.type == "thetas") {
    return(list(
      fixed = object$fixed.form,
      random = object$random.mu,
      start = object$theta
    ))
  } else if (mu.type == "covariates") {
    theta <- names(object$theta)
    cov.ref <- object$cov.ref
    cov.theta <- unique(as.vector(unlist(cov.ref)))
    cov.base <- theta[!(theta %in% cov.theta)]
    cov.base <- cov.base[!(cov.base %in% unlist(lapply(names(cov.ref), function(x) {
      names(cov.ref[[x]])
    })))]
    cov.lst <- list()
    new.theta <- cov.base
    for (n in names(cov.ref)) {
      cov.base <- cov.base[!(cov.base %in% (names(cov.ref[[n]])))]
      cur <- cov.ref[[n]]
      for (i in seq_along(cur)) {
        m <- cur[i]
        cov.lst[[m]] <- c(cov.lst[[m]], n)
        new.theta <- c(new.theta, as.vector(m), names(m))
      }
    }
    e1 <- paste(paste(cov.base, collapse = "+"), "~ 1")
    fixed.form <- paste(c(e1, sapply(names(cov.lst), function(x) {
      paste(x, "~", paste(cov.lst[[x]], collapse = "+"))
    })), collapse = ", ")
    fixed.form <- eval(parse(text = sprintf("list(%s)", fixed.form)))
    if (length(cov.base) == 0) {
      fixed.form <- fixed.form[-1]
    }
    theta <- theta[new.theta]
    return(list(
      fixed = fixed.form,
      random = object$random.mu,
      start = object$theta
    ))
  } else {
    return(list(
      fixed = object$fixed.form,
      random = object$random,
      start = object$theta
    ))
  }
}
##' Create the nlme parameter transform function from the UI object.
##'
##' @param object UI object
##' @param mu Is the model mu referenced?
##' \itemize{
##'
##' \item With the "thetas" only the population parameters are
##' mu-referenced; All covariates are included in the model parameter
##' function.  The between subject variability pieces are specified in
##' the \code{random} specs parameter.
##'
##' \item With the "covariates" option, the population parameters are
##' mu referenced and covariates are removed from the model function.
##' The covariates will be specified used in the fixed effects
##' parameterization of nlme, like \code{list(lKA+lCL~1, lV~WT)}
##'
##' \item With the "none" option, the model function is given to nlme
##' without any modification.
##'
##' }
##' @return Parameter function for nlme
##' @author Matthew L. Fidler
##' @keywords internal
nlmixr2UI.nlmefun <- function(object, mu.type = c("thetas", "covariates", "none")) {
  ## create nlme function
  mu.type <- match.arg(mu.type)
  if (mu.type == "thetas") {
    if (length(object$mu.ref) == 0L) {
      return(NULL)
    }
    .all <- object$ini$name[which(object$ini$neta1 == object$ini$neta2)] %in% names(object$mu.ref)
    if (!all(.all)) {
      return(NULL)
    }
    fn <- eval(parse(text = sprintf("function(%s) NULL", paste(unique(c(names(object$ini$theta), object$all.covs)), collapse = ", "))))
    body(fn) <- body(object$nlme.mu.fun)
  } else if (mu.type == "covariates") {
    if (length(object$mu.ref) == 0L) {
      return(NULL)
    }
    vars <- unique(c(unlist(object$mu.ref), unlist(object$cov.ref)))
    fn <- eval(parse(text = sprintf("function(%s) NULL", paste(vars, collapse = ", "))))
    body(fn) <- body(object$nlme.mu.fun2)
    vars2 <- allVars(body(fn))
    if (length(vars) != length(vars2)) {
      return(NULL)
    }
  } else {
    fn <- eval(parse(text = sprintf("function(%s) NULL", paste(object$rest.vars, collapse = ", "))))
    body(fn) <- body(object$rest)
  }
  return(fn)
}
##' Return dynmodel variable translation function
##'
##' @param object nlmixr2 ui object
##' @return nlmixr2 dynmodel translation
##' @author Matthew Fidler
nlmixr2UI.dynmodelfun <- function(object) {
  .fn <- nlmixr2UI.nlmefun(object, "none")
  .fn <- deparse(body(.fn))
  .fn[1] <- paste0("{\n.env <-environment();\nsapply(names(..par),function(x){assign(x,setNames(..par[x],NULL),envir=.env)})\n")
  .fn[length(.fn)] <- paste("return(unlist(as.list(environment())))}")
  .fn <- eval(parse(text = paste0("function(..par)", paste(.fn, collapse = "\n"))))
  return(.fn)
}

##' Return dynmodel variable translation function
##'
##' @param object nlmixr2 ui object
##' @return nlmixr2 dynmodel translation
##' @author Matthew Fidler
nlmixr2UI.dynmodelfun2 <- function(object) {
  .fn <- nlmixr2UI.nlmefun(object, "none")
  .bfn <- body(.fn)
  .fn <- deparse(.bfn)
  .extra <- paste(deparse(nlmixr2findLhs(.bfn)), collapse = " ")
  .fn[1] <- paste0("{\n.env <-environment();\nsapply(names(..par),function(x){assign(x,setNames(..par[[x]],NULL),envir=.env)})\n")
  .fn[length(.fn)] <- paste0(
    ".names <- unique(c(names(..par),",
    .extra, "));\nreturn(as.data.frame(setNames(lapply(.names,function(x){get(x,setNames(..par[x],NULL), envir=.env)}),.names)));\n}"
  )
  .fn <- eval(parse(text = paste0("function(..par)", paste(.fn, collapse = "\n"))))
  return(.fn)
}


##' Get the variance for the nlme fit process based on UI
##'
##' @param object UI object
##' @return nlme/lme variance object
##' @author Matthew L. Fidler
##' @keywords internal
nlmixr2UI.nlme.var <- function(object) {
  ## Get the variance for the nlme object
  add.prop.errs <- object$add.prop.errs
  w.no.add <- which(!add.prop.errs$add)
  w.no.prop <- which(!add.prop.errs$prop)
  const <- grp <- ""
  power <- ", fixed=c(1)"
  powera <- ", fixed=list(power=1)"
  if (length(add.prop.errs$y) > 1) {
    grp <- " | nlmixr2.grp"
  }
  if (length(w.no.add) > 0) {
    const <- sprintf(", fixed=list(%s)", paste(paste0(add.prop.errs$y[w.no.add], "=0"), collapse = ", "))
  }
  if (length(w.no.prop) > 0) {
    power <- sprintf(", fixed=list(%s)", paste(paste0(add.prop.errs$y, "=", ifelse(add.prop.errs$prop, 1, 0)), collapse = ", "))
    powera <- sprintf(", fixed=list(power=list(%s))", paste(paste0(add.prop.errs$y, "=", ifelse(add.prop.errs$prop, 1, 0)), collapse = ", "))
  }
  ## nlme_3.1-149 doesn't support varConstProp need to see if it exists
  .nlme <- loadNamespace("nlme")
  .varConstProp <- FALSE
  if (length(ls(pattern = "^varConstProp$", envir = .nlme)) == 1L) {
    .varConstProp <- TRUE
  }
  .addProp <- object$env$.addProp
  if (.addProp == "combined2" && !.varConstProp) {
    warning("this version of nlme does not support combined2 add+prop, degrading to combined1")
    object$env$.addProp <- .addProp <- "combined1"
  }
  if (.addProp == "combined1") {
    tmp <- sprintf("varConstPower(form=~fitted(.)%s%s)", grp, powera)
  } else if (.addProp == "combined2") {
    if (powera != ", fixed=list(power=1)") {
      stop("add+prop combined2 does not support nlme power currently")
    }
    tmp <- sprintf("nlme::varConstProp(%s)", grp)
  }
  if (all(!add.prop.errs$prop)) {
    tmp <- sprintf("varIdent(form = ~ 1%s)", grp)
    if (tmp == "varIdent(form = ~ 1)") {
      warning("initial condition for additive error ignored with nlme")
      return(NULL)
    }
  } else if (all(!add.prop.errs$add)) {
    tmp <- sprintf("varPower(form = ~ fitted(.)%s%s)", grp, power)
  }
  return(eval(parse(text = tmp)))
}
##' Return rxode2 model with predictions appended
##'
##' @param object UI object
##' @return String or NULL if rxode2 is not specified by UI.
##' @author Matthew L. Fidler
nlmixr2UI.rxode.pred <- function(object) {
  if (is.null(object$rxode)) {
    return(NULL)
  } else {
    tmp <- .deparse1(body(object$pred))
    return(paste(c(object$rxode, tmp), collapse = "\n"))
  }
}

##' Return rxode2 model with predictions appended
##'
##' @param object UI object
##' @return Combined focei model text for rxode2
##' @author Matthew L. Fidler
nlmixr2UI.saem.rx1 <- function(object) {
  .prd <- .deparse1(body(object$pred))
  .w <- any(regexpr("\\bnlmixr2_lincmt_pred\\b", .prd, perl = TRUE) != -1)
  if (length(.w) == 1) {
    .prd[.w] <- paste0("nlmixr2_lincmt_pred <- linCmt()\n", .prd[.w])
  }
  if (exists(".sf1", object$env)) {
    .bf <- .deparse1(body(object$env$.sf1))
  } else {
    .bf <- .deparse1(body(object$saem.fun1))
  }
  .ret <- paste(c(.bf, .prd, object$extra), collapse = "\n")
  if (any(regexpr("\\blinCmt[(]", .prd, perl = TRUE) != -1)) {
    .ret <- rxode2::rxNorm(rxode2::rxGetLin(.ret))
  }
  .ret
}

##' Return rxode2 model with predictions appended
##'
##' @param object UI object
##' @return Combined focei model text for rxode2
##' @author Matthew L. Fidler
##' @noRd
nlmixr2UI.focei.rx1 <- function(obj) {
  .df <- .as.data.frame(obj$ini)
  .dft <- .df[!is.na(.df$ntheta), ]
  .unfixed <- with(.dft, sprintf("%s=THETA[%d]", name, seq_along(.dft$name)))
  .eta <- .df[!is.na(.df$neta1), ]
  .eta <- .eta[.eta$neta1 == .eta$neta2, ]
  .eta <- with(.eta, sprintf("%s=ETA[%d]", name, .eta$neta1))
  .prd <- .deparse1(body(obj$pred))
  .w <- any(regexpr("\\bnlmixr2_lincmt_pred\\b", .prd, perl = TRUE) != -1)
  if (length(.w) == 1) {
    .prd[.w] <- paste0("nlmixr2_lincmt_pred <- linCmt()\n", .prd[.w])
  }
  .ret <- paste(c(
    .unfixed, .eta, .deparse1(body(obj$focei.fun1)),
    .prd, obj$extra
  ), collapse = "\n")
  if (any(regexpr("\\blinCmt[(]", .prd, perl = TRUE) != -1)) {
    .ret <- rxode2::rxNorm(rxode2::rxGetLin(.ret))
  }
  .ret
}

##' Get the Parameter  function with THETA/ETAs defined
##'
##' @param obj UI object
##' @return parameters function defined in THETA[#] and ETA[#]s.
##' @author Matthew L. Fidler
nlmixr2UI.theta.pars <- function(obj) {
  .df <- .as.data.frame(obj$ini)
  .dft <- .df[!is.na(.df$ntheta), ]
  .unfixed <- with(.dft, sprintf("%s=THETA[%d]", name, seq_along(.dft$name)))
  .eta <- .df[!is.na(.df$neta1), ]
  .eta <- .eta[.eta$neta1 == .eta$neta2, ]
  .eta <- with(.eta, sprintf("%s=ETA[%d]", name, .eta$neta1))
  .f <- .deparse1(body(obj$rest))
  .f <- eval(parse(text = paste(c("function(){", .unfixed, .eta, .f, "}"), collapse = "\n")))
  return(.f)
}
##' Get SAEM distribution
##'
##' @param obj UI object
##' @return Character of distribution
##' @author Matthew L. Fidler
nlmixr2UI.saem.distribution <- function(obj) {
  .df <- obj$ini$err
  .df <- paste(.df[which(!is.na(.df))])
  if (any(.df %in% c("dpois", "pois"))) {
    return("poisson")
  }
  if (any(.df %in% c("dbern", "bern", "dbinom", "binom"))) {
    if (.df %in% c("dbinom", "binom")) {
      .df <- obj$ini
      .w <- which(.df$err %in% c("dbinom", "binom"))
      if (length(.w) != 1L) stop("Distribution unsupported by SAEM")
      if (!is.na(.df$name[.w])) stop("Distribution unsupported by SAEM")
    }
    return("binomial")
  }
  if (any(.df %in% c("dlnorm", "lnorm", "logn", "dlogn", "logitNorm", "probitNorm"))) {
    return("normal")
  }
  if (any(.df %in% c("dnorm", "norm", "prop", "propT", "add", "pow", "powT", "pow2", "powT2"))) {
    return("normal")
  }

  stop("Distribution unsupported by SAEM")
}
##' Get parameters that are fixed
##'
##' @param obj UI object
##' @return logical vector of fixed THETA parameters
##' @author Matthew L. Fidler
nlmixr2UI.focei.fixed <- function(obj) {
  .df <- .as.data.frame(obj$ini)
  .dft <- .df[!is.na(.df$ntheta), ]
  .fix <- .dft$fix
  .dft <- .df[is.na(.df$ntheta), ]
  .fix <- c(.fix, .dft$fix)
  return(.fix)
}
##' Get parameters that are fixed for SAEM
##'
##' @param obj UI object
##' @return List of parameters that are fixed.
##' @author Matthew L. Fidler
nlmixr2UI.saem.fixed <- function(obj) {
  .df <- .as.data.frame(obj$ini)
  .dft <- .df[!is.na(.df$ntheta), ]
  .fixError <- .dft[!is.na(.dft$err), ]
  if (any(.fixError$fix)) {
    stop("Residuals cannot be fixed in SAEM.")
  }
  .dft <- .dft[is.na(.dft$err), ]
  .dft <- setNames(.dft$fix, paste(.dft$name))
  .dft <- .dft[obj$saem.theta.name]
  return(setNames(which(.dft), NULL))
}

##' Get the FOCEi initializations
##'
##' @param obj UI object
##' @return list with FOCEi style initializations
##' @author Matthew L. Fidler
nlmixr2UI.focei.inits <- function(obj) {
  df <- .as.data.frame(obj$ini)
  dft <- df[!is.na(df$ntheta), ]
  eta <- df[!is.na(df$neta1), ]
  len <- length(eta$name)
  cur.lhs <- character()
  cur.rhs <- numeric()
  ome <- character()
  for (i in seq_along(eta$name)) {
    last.block <- FALSE
    if (i == len) {
      last.block <- TRUE
    } else if (eta$neta1[i + 1] == eta$neta2[i + 1]) {
      last.block <- TRUE
    }
    if (eta$neta1[i] == eta$neta2[i]) {
      cur.lhs <- c(cur.lhs, sprintf("ETA[%d]", eta$neta1[i]))
      cur.rhs <- c(cur.rhs, eta$est[i])
      if (last.block) {
        ome[length(ome) + 1] <- sprintf(
          "%s ~ %s", paste(cur.lhs, collapse = " + "),
          paste(.deparse(cur.rhs), collapse = " ")
        )
        cur.lhs <- character()
        cur.rhs <- numeric()
      }
    } else {
      cur.rhs <- c(cur.rhs, eta$est[i])
    }
  }
  ome <- eval(parse(text = sprintf("list(%s)", paste(ome, collapse = ","))))
  return(list(
    THTA = dft$est,
    OMGA = ome
  ))
}
##' Get the eta->eta.trans for SAEM
##'
##' @param obj ui object
##' @return list of eta to eta.trans
##' @author Matthew L. Fidler
nlmixr2UI.saem.eta.trans <- function(obj) {
  eta.names <- obj$eta.names
  theta.names <- obj$theta.names
  theta.trans <- .saemThetaTrans(obj)
  mu.ref <- obj$mu.ref
  trans <- rep(NA, length(eta.names))
  for (i in seq_along(eta.names)) {
    ref <- mu.ref[[eta.names[i]]]
    if (!is.null(ref)) {
      w <- which(ref == theta.names)
      if (length(w) == 1) {
        trans[i] <- theta.trans[w]
      }
    }
  }
  ## Take out error terms
  if (any(is.na(trans))) {
    stop("Could not figure out the mu-referencing for this model.")
  }
  return(trans)
}
##' Get the SAEM model Omega
##'
##' @param obj UI model
##' @return SAEM model$omega spec
##' @author Matthew L. Fidler
nlmixr2UI.saem.model.omega <- function(obj) {
  dm <- sum(!is.na(.saemThetaTrans(obj)))
  et <- obj$saem.eta.trans
  mat <- matrix(rep(0, dm * dm), dm)
  etd <- which(!is.na(obj$neta1))
  for (i in etd) {
    mat[et[obj$neta1[i]], et[obj$neta2[i]]] <- mat[et[obj$neta2[i]], et[obj$neta1[i]]] <- 1
  }
  return(mat)
}

nlmixr2UI.saem.low <- function(obj) {
  .predDf <- obj$predDf
  .ini <- .as.data.frame(obj$ini)
  .ini <- .ini[!is.na(.ini$err), ]
  return(sapply(.predDf$cond, function(x) {
    .tmp <- .ini[which(.ini$condition == x), ]
    .w <- which(.tmp$err == "logitNorm")
    if (length(.w) == 1L) {
      return(.tmp$trLow[.w])
    }
    .w <- which(.tmp$err == "probitNorm")
    if (length(.w) == 1L) {
      return(.tmp$trLow[.w])
    }
    return(-Inf)
  }))
}

nlmixr2UI.saem.hi <- function(obj) {
  .predDf <- obj$predDf
  .ini <- .as.data.frame(obj$ini)
  .ini <- .ini[!is.na(.ini$err), ]
  return(sapply(.predDf$cond, function(x) {
    .tmp <- .ini[which(.ini$condition == x), ]
    .w <- which(.tmp$err == "logitNorm")
    if (length(.w) == 1L) {
      return(.tmp$trHi[.w])
    }
    .w <- which(.tmp$err == "probitNorm")
    if (length(.w) == 1L) {
      return(.tmp$trHi[.w])
    }
    return(Inf)
  }))
}

nlmixr2UI.saem.propT <- function(obj) {
  if (any(obj$saem.distribution == c("poisson", "binomial"))) {
    return(0L)
  }
  .predDf <- obj$predDf
  .ini <- .as.data.frame(obj$ini)
  .ini <- .ini[!is.na(.ini$err), ]
  return(sapply(.predDf$cond, function(x) {
    .tmp <- .ini[which(.ini$condition == x), ]
    if (any(.tmp$err == "propT") | any(.tmp$err == "powT")) {
      return(1L)
    } else {
      return(0L)
    }
  }))
}

# Get the TBS yj for SAEM
nlmixr2UI.saem.yj <- function(obj) {
  ## yj is for the number of endpoints in the nlmir model
  ## yj = 7 Yeo Johnson and  probit()
  ## yj = 6 probit()
  ## yj = 5 Yeo Johnson and logit()
  ## yj = 4 logitNorm
  ## yj = 3 logNorm
  ## yj = 2 norm()
  ## yj = 0 boxCox + Norm
  ## yj = 1 yeoJohnson + Norm
  if (any(obj$saem.distribution == c("poisson", "binomial"))) {
    return(2L)
  }
  .predDf <- obj$predDf
  .ini <- .as.data.frame(obj$ini)
  .ini <- .ini[!is.na(.ini$err), ]
  return(sapply(.predDf$cond, function(x) {
    .tmp <- .ini[which(.ini$condition == x), ]
    .boxCox <- any(.tmp$err == "boxCox")
    .yeoJohnson <- any(.tmp$err == "yeoJohnson")
    .hasAdd0 <- any(.tmp$err == "add") | any(.tmp$err == "norm") | any(.tmp$err == "dnorm")
    .hasLog <- any(.tmp$err == "dlnorm") | any(.tmp$err == "lnorm") | any(.tmp$err == "logn") |
      any(.tmp$err == "dlogn")
    .hasLogit <- any(.tmp$err == "logitNorm")
    .hasProbit <- any(.tmp$err == "probitNorm")
    if (.hasLog) {
      return(3L)
    }
    if (.hasLogit & .yeoJohnson) {
      return(5L)
    }
    if (.hasLogit) {
      return(4L)
    }
    if (.hasProbit & .yeoJohnson) {
      return(7L)
    }
    if (.hasProbit) {
      return(6L)
    }
    if (.boxCox) {
      return(0L)
    }
    if (.yeoJohnson) {
      return(1L)
    }
    return(2L)
  }))
}

# Get the lambda estimates for SAEM
nlmixr2UI.saem.lambda <- function(obj) {
  if (any(obj$saem.distribution == c("poisson", "binomial"))) {
    return(1.0)
  }
  .predDf <- obj$predDf
  .ini <- .as.data.frame(obj$ini)
  .ini <- .ini[!is.na(.ini$err), ]
  return(sapply(.predDf$cond, function(x) {
    .tmp <- .ini[which(.ini$condition == x), ]
    .boxCox <- which(.tmp$err == "boxCox")
    if (length(.boxCox) == 1L) {
      return(.tmp$est[.boxCox])
    }
    .yeoJohnson <- which(.tmp$err == "yeoJohnson")
    if (length(.yeoJohnson) == 1L) {
      return(.tmp$est[.yeoJohnson])
    }
    return(1.0)
  }))
}
##' Get the SAEM model$res.mod code
##'
##' @param obj UI model
##' @return SAEM model$res.mod spec
##' @author Matthew L. Fidler
nlmixr2UI.saem.res.mod <- function(obj) {
  if (any(obj$saem.distribution == c("poisson", "binomial"))) {
    return(1)
  }
  .predDf <- obj$predDf
  .ini <- .as.data.frame(obj$ini)
  .ini <- .ini[!is.na(.ini$err), ]
  return(sapply(.predDf$cond, function(x) {
    .tmp <- .ini[which(.ini$condition == x), ]
    .hasAdd0 <- any(.tmp$err == "add") | any(.tmp$err == "norm") | any(.tmp$err == "dnorm")
    .hasLog <- FALSE
    .w <- c(which(.tmp$err == "dlnorm"), which(.tmp$err == "lnorm"), which(.tmp$err == "logn"), which(.tmp$err == "dlogn"))
    if (length(.w) == 1) {
      if (!is.na(.tmp$est[.w])) {
        .hasLog <- TRUE
      }
    }
    .w <- which(.tmp$err == "logitNorm")
    .hasLogit <- FALSE
    if (length(.w) == 1) {
      if (!is.na(.tmp$est[.w])) {
        .hasLogit <- TRUE
      }
    }
    .hasProbit <- FALSE
    .w <- which(.tmp$err == "probitNorm")
    if (length(.w) == 1) {
      if (!is.na(.tmp$est[.w])) {
        .hasProbit <- TRUE
      }
    }
    .hasAdd <- .hasAdd0 | .hasLog | .hasLogit | .hasProbit
    .hasProp <- any(.tmp$err == "prop") | any(.tmp$err == "propT")
    .hasPow <- any(.tmp$err == "pow") | any(.tmp$err == "powT")
    .boxCox <- which(.tmp$err == "boxCox")
    .hasLambda <- FALSE
    if (length(.boxCox) == 1L) {
      if (!.tmp$fix[.boxCox]) .hasLambda <- TRUE
    }
    .yeoJohnson <- which(.tmp$err == "yeoJohnson")
    if (length(.yeoJohnson) == 1L) {
      if (!.tmp$fix[.yeoJohnson]) {
        if (.hasLambda) stop("cannot use both Yeo-Johnson and Box-Cox transformations", call. = FALSE)
        .hasLambda <- TRUE
      }
    }
    ## mod$res.mod = 10 = additive + pow + lambda
    if (.hasPow & .hasLambda & .hasAdd & !.hasProp) {
      return(10L)
    }
    ## mod$res.mod = 9 = additive + prop + lambda
    if (.hasProp & .hasLambda & .hasAdd & !.hasPow) {
      return(9L)
    }
    ## mod$res.mod = 8 = power + lambda
    if (.hasPow & .hasLambda & !.hasAdd & !.hasProp) {
      return(8L)
    }
    ## mod$res.mod = 7 = prop + lambda
    if (.hasProp & .hasLambda & !.hasAdd & !.hasPow) {
      return(7L)
    }
    ## mod$res.mod = 6 = additive + lambda
    if (.hasAdd & .hasLambda & !.hasProp & !.hasPow) {
      return(6L)
    }
    ## mod$res.mod = 5 = power
    if (.hasPow & !.hasAdd & !.hasProp & !.hasLambda) {
      return(5L)
    }
    ## mod$res.mod = 4 = additive + power
    if (.hasAdd & .hasPow & !.hasLambda & !.hasProp) {
      return(4L)
    }
    ## mod$res.mod = 3 = additive + proportional
    if (.hasAdd & .hasProp & !.hasLambda & !.hasPow) {
      return(3L)
    }
    ## mod$res.mod = 2 = proportional
    if (.hasProp & !.hasPow & !.hasAdd & !.hasLambda) {
      return(2L)
    }
    ## mod$res.mod = 1 = additive or poisson
    if (.hasAdd & !.hasPow & !.hasProp & !.hasLambda) {
      return(1L)
    }
  }))
}
##' Get error names for SAEM
##'
##' @param obj SAEM user interface function.
##' @return Names of error estimates for SAEM
##' @author Matthew L. Fidler
nlmixr2UI.saem.res.name <- function(obj) {
  w <- which(sapply(obj$err, function(x) any(x == c("add", "norm", "dnorm", "dlnorm", "lnorm", "logn", "dlogn"))))
  ret <- c()
  if (length(w) == 1) {
    if (!is.na(obj$est[w])) {
      ret[length(ret) + 1] <- paste(obj$name[w])
    }
  }
  w <- c(which(obj$err == "prop"), which(obj$err == "propT"))
  if (length(w) == 1) {
    ret[length(ret) + 1] <- paste(obj$name[w])
  }
  return(ret)
}

##' Get initial estimate for ares SAEM.
##'
##' @param obj UI model
##' @return SAEM model$ares spec
##' @author Matthew L. Fidler
nlmixr2UI.saem.ares <- function(obj) {
  .predDf <- obj$predDf
  .ini <- .as.data.frame(obj$ini)
  .ini <- .ini[!is.na(.ini$err), ]
  return(sapply(.predDf$cond, function(x) {
    .tmp <- .ini[which(.ini$condition == x), ]
    .w <- which(sapply(.tmp$err, function(x) {
      any(x == c(
        "add", "norm", "dnorm", "dpois",
        "pois", "dbinom", "binom", "dbern", "bern",
        "lnorm", "dlnorm", "logn", "dlogn"
      ))
    }))
    if (length(.w) == 1) {
      return(.tmp$est[.w])
    } else {
      return(10)
    }
  }))
}

##' Get initial estimate for bres SAEM.
##'
##' @param obj UI model
##' @return SAEM model$ares spec
##' @author Matthew L. Fidler
nlmixr2UI.saem.bres <- function(obj) {
  .predDf <- obj$predDf
  .ini <- .as.data.frame(obj$ini)
  .ini <- .ini[!is.na(.ini$err), ]
  return(sapply(.predDf$cond, function(x) {
    .tmp <- .ini[which(.ini$condition == x), ]
    .w <- which(sapply(.tmp$err, function(x) (any(x == "prop") || any(x == "propT"))))
    if (length(.w) == 1) {
      return(.tmp$est[.w])
    } else {
      .w <- which(sapply(.tmp$err, function(x) (any(x == "pow") || any(x == "powT"))))
      if (length(.w) == 1) {
        return(.tmp$est[.w])
      } else {
        return(1)
      }
    }
  }))
}

##' Get initial estimate for bres SAEM.
##'
##' @param obj UI model
##' @return SAEM model$ares spec
##' @author Matthew L. Fidler
nlmixr2UI.saem.cres <- function(obj) {
  .predDf <- obj$predDf
  .ini <- .as.data.frame(obj$ini)
  .ini <- .ini[!is.na(.ini$err), ]
  return(sapply(.predDf$cond, function(x) {
    .tmp <- .ini[which(.ini$condition == x), ]
    .w <- which(sapply(.tmp$err, function(x) (any(x == "pow2") || any(x == "powT2"))))
    if (length(.w) == 1) {
      return(.tmp$est[.w])
    } else {
      return(1)
    }
  }))
}

##' Get model$log.eta for SAEM
##'
##' @param obj UI model
##' @return SAEM model$log.eta
##' @author Matthew L. Fidler
nlmixr2UI.saem.log.eta <- function(obj) {
  lt <- obj$log.theta
  dm <- sum(!is.na(.saemThetaTrans(obj)))
  ret <- rep(FALSE, dm)
  theta.trans <- .saemThetaTrans(obj)
  theta.names <- obj$theta.names
  for (n in lt) {
    w <- which(n == theta.names)
    if (length(w) == 1) {
      ret[theta.trans[w]] <- TRUE
    }
  }
  return(ret)
}

##' Generate saem.fit user function.
##'
##' @param obj UI object
##' @return saem user function
##' @author Matthew L. Fidler
nlmixr2UI.saem.fit <- function(obj) {
  if (any(ls(envir = obj$env) == "saem.fit")) {
    return(obj$env$saem.fit)
  } else if (!is.null(obj$rxode)) {
    ## rxode2 function
    if (is.null(obj$saem.pars)) {
      stop("SAEM requires mu-referenced parameters")
    }
    inPars <- obj$saem.inPars
    if (length(inPars) == 0) {
      inPars <- NULL
    } else {
      ## Check for inPars in Covariates in rxode2 model
      .extra <- rxode2::rxModelVars(obj$rxode)
      .extra <- intersect(obj$saem.all.covs, .extra$params)
      if (length(.extra) > 0) {
        inPars <- c(inPars, .extra)
      }
    }
    pars <- obj$saem.pars
    if (!is.null(obj$env$.curTv)) {
      ## need to add back .curTv to pars
      .curTv <- obj$env$.curTv
      .covRef <- obj$nmodel$cov.ref
      .f <- function(x, theta, extra) {
        if (is.atomic(x)) {
          return(x)
        } else if (is.name(x)) {
          if (any(as.character(x) == theta)) {
            return(eval(parse(
              text = paste0("quote(", x, "+", extra, ")")
            )))
          } else {
            return(x)
          }
        } else if (is.call(x)) {
          return(as.call(lapply(x, .f, theta = theta, extra = extra)))
        }
      }
      .bpars <- body(pars)
      .sf1 <- obj$saem.fun1
      .bfun1 <- body(.sf1)
      .extra <- c()
      for (.var in .curTv) {
        .lst <- .covRef[[.var]]
        .extra <- c(.extra, names(.lst))
        .estPar <- paste0(.var, "*", names(.lst))
        .thetaPar <- setNames(.lst, NULL)
        for (i in seq_along(.estPar)) {
          .bpars <- .f(.bpars, .thetaPar[i], .estPar[i])
          .bfun1 <- .f(.bfun1, .thetaPar[i], .estPar[i])
        }
      }
      body(pars) <- .bpars
      body(.sf1) <- .bfun1
      obj$env$.bpars <- pars
      inPars <- unique(c(inPars, .curTv))
      ## Also modify saem.rx1
      obj$env$.sf1 <- .sf1
    }
    ## saem.fit <- gen_saem_user_fn(model = ode, obj$saem.pars, pred = obj$predSaem, inPars = inPars)
    if (!exists("singleOde", obj$env)) {
      obj$env$singleOde <- TRUE
    }
    if (obj$env$singleOde) {
      saem.fit <- gen_saem_user_fn(obj$saem.rx1, NULL, obj$predSaem, inPars = inPars)
    } else {
      saem.fit <- gen_saem_user_fn(obj$rxode.pred, pars, obj$predSaem, inPars = inPars)
    }
    ## obj$env$saem.ode <- attr(saem.fit, "rx")
    obj$env$saem.fit <- saem.fit
    return(obj$env$saem.fit)
  }
}
##' Generate SAEM model list
##'
##' @param obj  nlmixr2 UI object
##' @return SAEM model list
##' @author Matthew L. Fidler
nlmixr2UI.saem.model <- function(obj) {
  mod <- list(saem_mod = obj$saem.fit)
  if (length(obj$saem.all.covs > 0)) {
    if (!is.null(obj$env$.curTv)) {
      .covars <- setdiff(obj$saem.all.covs, obj$env$.curTv)
      if (length(.covars) > 0) {
        mod$covars <- .covars
      }
    } else {
      mod$covars <- obj$saem.all.covs
    }
  }
  mod$res.mod <- obj$saem.res.mod
  mod$log.eta <- obj$saem.log.eta
  ## if (FALSE){
  ## FIXME option/warning
  mod$ares <- obj$saem.ares
  mod$bres <- obj$saem.bres
  ## return(nlmixr2UI.saem.bres(obj))
  ## }
  mod$omega <- obj$saem.model.omega
  return(mod)
}
.saemThetaTrans <- function(uif, muOnly = FALSE) {
  .tv <- uif$env$.curTv
  if (is.null(.tv)) {
    .trans <- uif$saem.theta.trans
  } else {
    .ret <- as.character(body(uif$env$.bpars))
    if (.ret[1] == "{") .ret <- .ret[-1]
    .pars <- rxode2::rxModelVars(paste(.ret, collapse = "\n"))$params
    .thetaNames <- uif$ini$theta.names
    .pars <- setdiff(.pars, uif$all.covs)
    .trans <- setNames(sapply(.thetaNames, function(x) {
      .w <- which(x == .pars)
      if (length(.w) == 1) {
        return(.w)
      }
      return(NA_integer_)
    }), NULL)
  }
  if (muOnly) {
    .trans <- setdiff(.trans, uif$ini$ntheta[which(!is.na(uif$ini$err) & !is.na(uif$ini$ntheta))])
  }
  return(.trans)
}
.saemAllCov <- function(uif) {
  all.covs <- uif$saem.all.covs
  .tv <- uif$env$.curTv
  if (!is.null(.tv)) {
    all.covs <- setdiff(all.covs, uif$env$.curTv)
  }
  return(all.covs)
}
##' Get THETA names for nlmixr2's SAEM
##'
##' @param uif nlmixr2 UI object
##' @return SAEM theta names
##' @author Matthew L. Fidler
nlmixr2UI.saem.theta.name <- function(uif) {
  .trans <- .saemThetaTrans(uif)
  .df <- .as.data.frame(uif$ini)
  .df <- .df[!is.na(.df$ntheta), ]
  .transName <- paste(.df$name[which(!is.na(.trans))])
  .trans <- .trans[!is.na(.trans)]
  theta.name <- .transName[order(.trans)]
  all.covs <- .saemAllCov(uif)
  lc <- length(all.covs)
  if (lc > 0) {
    m <- matrix(rep(NA, length(theta.name) * (lc + 1)), nrow = lc + 1)
    dimnames(m) <- list(c("_name", all.covs), theta.name)
    m["_name", ] <- theta.name
    for (cn in names(uif$cov.ref)) {
      v <- uif$cov.ref[[cn]]
      for (var in names(v)) {
        rn <- v[var]
        m[cn, rn] <- var
      }
    }
    ret <- unlist(m)
    ret <- ret[!is.na(ret)]
    return(ret)
  }
  return(theta.name)
}
##' Generate SAEM initial estimates for THETA.
##'
##' @param obj nlmixr2 UI object
##' @return SAEM theta initial estimates
##' @author Matthew L. Fidler
nlmixr2UI.saem.init.theta <- function(obj) {
  theta.name <- obj$saem.theta.name
  .tv <- obj$env$.curTv
  if (is.null(.tv)) {
    cov.names <- unique(names(unlist(structure(obj$cov.ref, .Names = NULL))))
  } else {
    cov.names <- unique(names(unlist(structure(obj$cov.ref[!(names(obj$cov.ref) %in% .tv)], .Names = NULL))))
  }
  theta.name <- theta.name[!(theta.name %in% cov.names)]
  nm <- paste(obj$ini$name)
  lt <- obj$log.theta
  i <- 0
  this.env <- environment()
  theta.ini <- sapply(theta.name, function(x) {
    w <- which(x == nm)
    assign("i", i + 1, this.env)
    if (any(lt == x)) {
      return(exp(obj$ini$est[w]))
    } else {
      return(obj$ini$est[w])
    }
  })
  all.covs <- .saemAllCov(obj)
  lc <- length(all.covs)
  if (lc > 0) {
    m <- matrix(rep(NA, lc * length(theta.name)), ncol = lc)
    dimnames(m) <- list(theta.name, all.covs)
    for (cn in names(obj$cov.ref)) {
      if (!(cn %in% .tv)) {
        v <- obj$cov.ref[[cn]]
        for (var in names(v)) {
          rn <- v[var]
          w <- which(var == nm)
          m[rn, cn] <- obj$ini$est[w]
        }
      }
    }
    return(as.vector(c(theta.ini, as.vector(m))))
  }
  return(as.vector(theta.ini))
}
##' SAEM's init$omega
##'
##' @param obj nlmixr2 UI object
##' @param names When \code{TRUE} return the omega names.  By default
##'     this is \code{FALSE}.
##' @return Return initial matrix
##' @author Matthew L. Fidler
nlmixr2UI.saem.init.omega <- function(obj, names = FALSE) {
  dm <- sum(!is.na(.saemThetaTrans(obj)))
  et <- obj$saem.eta.trans
  ret <- rep(NA, dm)
  etd <- which(obj$neta1 == obj$neta2)
  for (i in etd) {
    if (names) {
      ret[et[obj$neta1[i]]] <- paste(obj$name[i])
    } else {
      ret[et[obj$neta1[i]]] <- obj$est[i]
    }
  }
  if (names) {
    ret <- ret[!is.na(ret)]
    return(ret)
  } else {
    tmp <- unique(ret[!is.na(ret)])
    if (length(tmp) == 1) {
      ret[is.na(ret)] <- tmp
    } else {
      ret[is.na(ret)] <- 1
    }
  }
  return(ret)
}
##' Get saem initilization list
##'
##' @param obj nlmixr2 UI object
##' @return Return SAEM inits list.
##' @author Matthew L. Fidler
nlmixr2UI.saem.init <- function(obj) {
  ret <- list()
  ret$theta <- obj$saem.init.theta
  ## if (FALSE){
  ret$omega <- obj$saem.init.omega
  ## }
  return(ret)
}

nlmixr2UI.focei.mu.ref <- function(obj) {
  .muRef <- obj$mu.ref
  .tn <- obj$focei.names
  sapply(seq_along(.muRef), function(x) {
    .cur <- .muRef[[x]]
    .w <- which(.tn == .cur)
    if (length(.w) == 1) {
      return(.w - 1)
    } else {
      return(-1)
    }
  })
}

nlmixr2UI.model.desc <- function(obj) {
  .mv <- rxode2::rxModelVars(obj$rxode.pred)
  if (obj$predSys) {
    return("rxode2-based Pred model")
  } else if (.mv$extraCmt == 0) {
    return("rxode2-based ODE model")
  } else {
    return(sprintf(
      "rxode2-based %s-compartment model%s", .mv$flags["ncmt"],
      ifelse(.mv$extraCmt == 2, " with first-order absorption", "")
    ))
  }
}

nlmixr2UI.poped.notfixed_bpop <- function(obj) {
  .df <- .as.data.frame(obj$ini)
  .tmp <- .df[!is.na(.df$ntheta) & is.na(.df$err), ]
  return(setNames(1 - .tmp$fix * 1, paste(.tmp$name)))
}

nlmixr2UI.poped.d <- function(obj) {
  .df <- .as.data.frame(obj$ini)
  .tmp <- .df[which(is.na(.df$ntheta) & .df$neta1 == .df$neta2), ]
  return(setNames(.tmp$est, paste(.tmp$name)))
}

nlmixr2UI.poped.sigma <- function(obj) {
  .df <- .as.data.frame(obj$ini)
  .tmp <- .df[!which(is.na(.df$err) & .df$neta1 == .df$neta2), ]
  return(setNames(.tmp$est * .tmp$est, paste(.tmp$name)))
}

nlmixr2UI.logThetasList <- function(obj) {
  .ini <- .as.data.frame(obj$ini)
  .logThetas <- as.integer(which(setNames(sapply(obj$focei.names, function(x) any(x == obj$log.theta)), NULL)))
  .thetas <- .ini[!is.na(.ini$ntheta), ]
  .one <- obj$oneTheta
  .logThetasF <- .thetas[.thetas$name %in% .one, "ntheta"]
  .logThetasF <- intersect(.logThetas, .logThetasF)
  list(.logThetas, .logThetasF)
}

nlmixr2UI.logitThetasList <- function(obj) {
  .ini <- .as.data.frame(obj$ini)
  .logitThetas <- as.integer(which(setNames(sapply(obj$focei.names,
                                                   function(x) any(x == obj$logit.theta)), NULL)))
  .thetas <- .ini[!is.na(.ini$ntheta), ]
  .one <- obj$oneTheta
  .logitThetasF <- .thetas[.thetas$name %in% .one, "ntheta"]
  .logitThetasF <- intersect(.logitThetas, .logitThetasF)
  list(.logitThetas, .logitThetasF)
}

nlmixr2UI.logitThetasListHi <- function(obj) {
  .thetaList <- nlmixr2UI.logitThetasList(obj)
  .names <- obj$focei.names
  .hi <- obj$logit.theta.hi
  list(
    .hi[.names[.thetaList[[1]]]],
    .hi[.names[.thetaList[[2]]]]
  )
}
nlmixr2UI.logitThetasListLow <- function(obj) {
  .thetaList <- nlmixr2UI.logitThetasList(obj)
  .names <- obj$focei.names
  .low <- obj$logit.theta.low
  list(
    .low[.names[.thetaList[[1]]]],
    .low[.names[.thetaList[[2]]]]
  )
}


nlmixr2UI.probitThetasList <- function(obj) {
  .ini <- .as.data.frame(obj$ini)
  .probitThetas <- as.integer(which(setNames(sapply(obj$focei.names, function(x) any(x == obj$probit.theta)), NULL)))
  .thetas <- .ini[!is.na(.ini$ntheta), ]
  .one <- obj$oneTheta
  .probitThetasF <- .thetas[.thetas$name %in% .one, "ntheta"]
  .probitThetasF <- intersect(.probitThetas, .probitThetasF)
  list(.probitThetas, .probitThetasF)
}

nlmixr2UI.probitThetasListHi <- function(obj) {
  .thetaList <- nlmixr2UI.probitThetasList(obj)
  .names <- obj$focei.names
  .hi <- obj$probit.theta.hi
  list(
    .hi[.names[.thetaList[[1]]]],
    .hi[.names[.thetaList[[2]]]]
  )
}
nlmixr2UI.probitThetasListLow <- function(obj) {
  .thetaList <- nlmixr2UI.probitThetasList(obj)
  .names <- obj$focei.names
  .low <- obj$probit.theta.low
  list(
    .low[.names[.thetaList[[1]]]],
    .low[.names[.thetaList[[2]]]]
  )
}

nlmixr2UI.poped.ff_fun <- function(obj) {
  if (!is.null(obj$lin.solved)) {
    stop("Solved system not supported yet.")
  } else {
    .df <- .as.data.frame(obj$ini)
    .dft <- .df[!is.na(.df$ntheta) & is.na(.df$err), ]
    .unfixed <- with(.dft, sprintf("%s=bpop[%d]", name, seq_along(.dft$name)))
    .eta <- .df[!is.na(.df$neta1), ]
    .eta <- .eta[.eta$neta1 == .eta$neta2, ]
    .eta <- with(.eta, sprintf("%s=b[%d]", name, .eta$neta1))
    .lhs <- nlmixr2findLhs(body(obj$rest))
    .f <- .deparse(body(obj$rest))[-1]
    .lhs <- sprintf("return(c(%s))", paste(sprintf("\"%s\"=%s", .lhs, .lhs), collapse = ", "))
    .f <- eval(parse(text = paste(c("function(x,a,bpop,b,bocc){", .unfixed, .eta, .f[-length(.f)], .lhs, "}"), collapse = "\n")))
    return(.f)
  }
}

##' @export
`$.nlmixr2UI` <- function(obj, arg, exact = TRUE) {
  x <- obj
  .cls <- class(x)
  class(x) <- "list"
  if (arg == "ini") {
    return(x$ini)
  } else if (arg == "nmodel") {
    return(x$nmodel)
  } else if (arg == "model") {
    return(x$model)
  } else if (arg == "nlme.fun.mu") {
    return(nlmixr2UI.nlmefun(obj, "thetas"))
  } else if (arg == "dynmodel.fun") {
    return(nlmixr2UI.dynmodelfun(obj))
  } else if (arg == "dynmodel.fun.df") {
    return(nlmixr2UI.dynmodelfun2(obj))
  } else if (arg == "nlme.fun") {
    return(nlmixr2UI.nlmefun(obj, "none"))
  } else if (arg == "nlme.fun.mu.cov") {
    return(nlmixr2UI.nlmefun(obj, "covariates"))
  } else if (arg == "nlme.specs") {
    return(nlmixr2UI.nlme.specs(obj, "none"))
  } else if (arg == "nlme.specs.mu") {
    return(nlmixr2UI.nlme.specs(obj, "thetas"))
  } else if (arg == "nlme.specs.mu.cov") {
    return(nlmixr2UI.nlme.specs(obj, "covariates"))
  } else if (arg == "nlme.var") {
    return(nlmixr2UI.nlme.var(obj))
  } else if (arg == "rxode.pred") {
    return(nlmixr2UI.rxode.pred(obj))
  } else if (arg == "theta.pars") {
    return(nlmixr2UI.theta.pars(obj))
  } else if (arg == "focei.inits") {
    return(nlmixr2UI.focei.inits(obj))
  } else if (arg == "focei.fixed") {
    return(nlmixr2UI.focei.fixed(obj))
  } else if (arg == "focei.mu.ref") {
    return(nlmixr2UI.focei.mu.ref(obj))
  } else if (arg == "saem.fixed") {
    return(nlmixr2UI.saem.fixed(obj))
  } else if (arg == "saem.theta.name") {
    return(nlmixr2UI.saem.theta.name(obj))
  } else if (arg == "saem.eta.trans") {
    return(nlmixr2UI.saem.eta.trans(obj))
  } else if (arg == "saem.model.omega") {
    return(nlmixr2UI.saem.model.omega(obj))
  } else if (arg == "saem.res.mod") {
    return(nlmixr2UI.saem.res.mod(obj))
  } else if (arg == "saem.ares") {
    return(nlmixr2UI.saem.ares(obj))
  } else if (arg == "saem.bres") {
    return(nlmixr2UI.saem.bres(obj))
  } else if (arg == "saem.cres") {
    return(nlmixr2UI.saem.cres(obj))
  } else if (arg == "saem.log.eta") {
    return(nlmixr2UI.saem.log.eta(obj))
  } else if (arg == "saem.yj") {
    return(nlmixr2UI.saem.yj(obj))
  } else if (arg == "saem.propT") {
    return(nlmixr2UI.saem.propT(obj))
  } else if (arg == "saem.hi") {
    return(nlmixr2UI.saem.hi(obj))
  } else if (arg == "saem.low") {
    return(nlmixr2UI.saem.low(obj))
  } else if (arg == "saem.lambda") {
    return(nlmixr2UI.saem.lambda(obj))
  } else if (arg == "saem.fit") {
    return(nlmixr2UI.saem.fit(obj))
  } else if (arg == "saem.model") {
    return(nlmixr2UI.saem.model(obj))
  } else if (arg == "saem.init.theta") {
    return(nlmixr2UI.saem.init.theta(obj))
  } else if (arg == "saem.init.omega") {
    return(nlmixr2UI.saem.init.omega(obj))
  } else if (arg == "saem.init") {
    return(nlmixr2UI.saem.init(obj))
  } else if (arg == "saem.omega.name") {
    return(nlmixr2UI.saem.init.omega(obj, TRUE))
  } else if (arg == "saem.res.name") {
    return(nlmixr2UI.saem.res.name(obj))
  } else if (arg == "model.desc") {
    return(nlmixr2UI.model.desc(obj))
  } else if (arg == "meta") {
    return(x$meta)
  } else if (arg == "saem.distribution") {
    return(nlmixr2UI.saem.distribution(obj))
  } else if (arg == "notfixed_bpop" || arg == "poped.notfixed_bpop") {
    return(nlmixr2UI.poped.notfixed_bpop(obj))
  } else if (arg == "poped.ff_fun") {
    return(nlmixr2UI.poped.ff_fun(obj))
  } else if (arg == "poped.d") {
    return(nlmixr2UI.poped.d(obj))
  } else if (arg == "poped.sigma") {
    return(nlmixr2UI.poped.sigma(obj))
  } else if (arg == "logThetasList") {
    return(nlmixr2UI.logThetasList(obj))
  } else if (arg == "logitThetasListLow") {
    return(nlmixr2UI.logitThetasListLow(obj))
  } else if (arg == "logitThetasListHi") {
    return(nlmixr2UI.logitThetasListHi(obj))
  } else if (arg == "logitThetasList") {
    return(nlmixr2UI.logitThetasList(obj))
  } else if (arg == "probitThetasListLow") {
    return(nlmixr2UI.probitThetasListLow(obj))
  } else if (arg == "probitThetasListHi") {
    return(nlmixr2UI.probitThetasListHi(obj))
  } else if (arg == "probitThetasList") {
    return(nlmixr2UI.probitThetasList(obj))
  } else if (arg == "focei.rx1") {
    return(nlmixr2UI.focei.rx1(obj))
  } else if (arg == "saem.rx1") {
    return(nlmixr2UI.saem.rx1(obj))
  } else if (arg == ".clean.dll") {
    if (exists(".clean.dll", envir = x$meta)) {
      clean <- x$meta$.clean.dll
      if (is(clean, "logical")) {
        return(clean)
      }
    }
    return(TRUE)
  } else if (arg == "random.mu") {
    return(nlmixr2BoundsOmega(x$ini, x$nmodel$mu.ref))
  } else if (arg == "bpop") {
    arg <- "theta"
  } else if (arg == "multipleEndpoint") {
    return(nlmixr2UI.multipleEndpoint(x))
  } else if (arg == "muRefTable") {
    class(x) <- .cls
    return(.nmMuTable(x))
  } else if (arg == "single.inner.1") {
    return(nlmixr2UI.inner.model(obj, TRUE, "combined1"))
  } else if (arg == "single.inner.2") {
    return(nlmixr2UI.inner.model(obj, TRUE, "combined2"))
  } else if (arg == "inner") {
    return(nlmixr2UI.inner.model(obj, TRUE, "combined2"))
  } else if (arg == "inner.par0") {
    return(nlmixr2UI.par0(obj))
  } else if (arg == "inner.numeric") {
    return(nlmixr2UI.inner.model(obj, TRUE, "combined2", only.numeric = TRUE))
  } else if (arg == "single.saem") {
    return(nlmixr2UI.saem.1(obj))
  } else if (arg == "pars.saem") {
    return(nlmixr2UI.saem.2(obj))
  } else if (arg == "single.saem.params"){
    return(nlmixr2UI.saem.1.params(obj))
  }
  m <- x$ini
  ret <- `$.nlmixr2Bounds`(m, arg, exact = exact)
  if (is.null(ret)) {
    m <- x$nmodel
    ret <- m[[arg, exact = exact]]
    if (is.null(ret)) {
      if (exists(arg, envir = x$meta)) {
        ret <- get(arg, envir = x$meta)
      }
    }
  }
  ret
}



##' @export
str.nlmixr2UI <- function(object, ...) {
  obj <- object
  class(obj) <- "list"
  str(obj$ini)
  str(obj$nmodel)
  cat(" $ ini       : Model initilizations/bounds object\n")
  cat(" $ model     : Original Model\n")
  cat(" $ model.desc: Model description\n")
  cat(" $ nmodel    : Parsed Model List\n")
  cat(" $ nlme.fun  : The nlme model function.\n")
  cat(" $ nlme.specs: The nlme model specs.\n")
  cat(" $ nlme.var  : The nlme model varaince.\n")
  cat(" $ rxode.pred: The rxode2 block with pred attached (final pred is nlmixr2_pred)\n")
  cat(" $ theta.pars: Parameters in terms of THETA[#] and ETA[#]\n")
  cat(" $ focei.inits: Initialization for FOCEi style blocks\n")
  cat(" $ focei.fixed: Logical vector of FOCEi fixed parameters\n")
  cat(" $ focei.mu.ref: Integer Vector of focei.mu.ref\n")
  cat(" $ saem.eta.trans: UI ETA -> SAEM ETA\n")
  cat(" $ saem.model.omega: model$omega for SAEM\n")
  cat(" $ saem.res.mod: model$res.mod for SAEM\n")
  cat(" $ saem.ares: model$ares for SAEM\n")
  cat(" $ saem.bres: model$bres for SAEM\n")
  cat(" $ saem.log.eta: model$log.eta for SAEM\n")
  cat(" $ saem.fit  : The SAEM fit user function\n")
  cat(" $ saem.model: The SAEM model list\n")
  cat(" $ saem.init.theta: The SAEM init$theta\n")
  cat(" $ saem.init.omega: The SAEM init$omega\n")
  cat(" $ saem.init : The SAEM inits list\n")
  cat(" $ saem.theta.name : The SAEM theta names\n")
  cat(" $ saem.omega.name : The SAEM theta names\n")
  cat(" $ saem.res.name : The SAEM omega names\n")
  cat(" $ saem.distribution: SAEM distribution\n")
  cat(" $ .clean.dll : boolean representing if dlls are cleaned after running.\n")
  cat(" $ logThetasList: List of logThetas:\n")
  cat("     first element are scaling log thetas;\n")
  cat("     second element are back-transformed thetas;\n")
  cat(" $ multipleEndpoint: table/huxtable of multiple endpoint translations in nlmixr2\n")
  cat(" $ muRefTable: table/huxtable of mu-referenced items in a model\n")
}

.syncUif <- function(uif, popDf = NULL, omega = NULL) {
  if (inherits(uif, "nlmixr2FitCore")) {
    popDf <- uif$popDf
    omega <- uif$omega
    uif <- uif$uif
  }
  .rn <- rownames(popDf)
  .est <- popDf$Estimate
  .ini <- as.data.frame(uif$ini)
  for (i in seq_along(popDf)) {
    .name <- .rn[i]
    .curEst <- .est[i]
    .w <- which(.ini$name == .name)
    .ini$est[.w] <- .curEst
  }
  .dn <- dimnames(omega)[[1]]
  for (.n1 in .dn) {
    .w <- which(.ini$name == .n1)
    .ini$est[.w] <- omega[.n1, .n1]
    for (.n2 in .dn) {
      .name <- paste0("(", .n1, ",", .n2, ")")
      .w <- which(.ini$name == .name)
      if (length(.w) == 1) {
        .ini$est[.w] <- omega[.n1, .n2]
      }
    }
  }
  class(.ini) <- c("nlmixr2Bounds", "data.frame")
  uif$ini <- .ini
  return(uif)
}

nlmixr2UI.inner.model <- function(obj, singleOde = TRUE, addProp = c("combined1", "combined2"), optExpression = TRUE,
                                 only.numeric = FALSE) {
  addProp <- match.arg(addProp)
  .mod <- obj$focei.rx1
  .pars <- NULL
  if (singleOde) {
    .mod <- obj$focei.rx1
    .pars <- NULL
  } else {
    .mod <- obj$rxode.pred
    .pars <- obj$theta.pars
  }
  if (missing(only.numeric)) {
    .onlyNumeric <- all(is.na(obj$ini$neta1))
  } else {
    .onlyNumeric <- only.numeric
  }
  rxode2::rxSymPySetupPred(.mod, function() {
    return(nlmixr2_pred)
  }, .pars, obj$error,
  grad = FALSE, pred.minus.dv = TRUE, sum.prod = FALSE,
  interaction = TRUE, only.numeric = .onlyNumeric, run.internal = TRUE,
  addProp = addProp, optExpression = optExpression
  )
}

nlmixr2UI.saem.1 <- function(obj, optExpression=FALSE, loadSymengine=FALSE) {
  .mod <- rxode2::rxode2(rxode2::rxGenSaem(obj$saem.rx1, function() {
    return(nlmixr2_pred)
  }, NULL,
  optExpression = optExpression,
  loadSymengine=loadSymengine))
  .mod
}

nlmixr2UI.saem.2 <- function(obj, optExpression=FALSE, loadSymengine=FALSE) {
  .mod <- rxode2::rxode2(rxode2::rxGenSaem(obj$rxode.pred, function() {
    return(nlmixr2_pred)
  }, obj$saem.pars, optExpression = optExpression, loadSymengine=loadSymengine))
  .mod
}

nlmixr2UI.saem.1.params <- function(obj) {
  .m <- nlmixr2UI.saem.1(obj)
  .df <- as.data.frame(obj$ini)
  .tmp <- sapply(.m$params, function(x){
    .w <- which(.df$name == x)
    if (length(.w) == 1) return(.df$est[.w])
    return(NA_real_)
  })
  .tmp <- .tmp[!is.na(.tmp)]
  .tmp
}

nlmixr2UI.par0 <- function(obj) {
  .ini <- as.data.frame(obj$ini)
  .maxEta <- max(.ini$neta1, na.rm = TRUE)
  .theta <- .ini[!is.na(.ini$ntheta), ]
  setNames(
    c(.theta$est, rep(0.0, .maxEta)),
    c(paste0("THETA[", .theta$ntheta, "]"), paste0("ETA[", seq(1, .maxEta), "]"))
  )
}


## Local Variables:
## ess-indent-offset: 2
## indent-tabs-mode: nil
## End:
