.fixNames <- function(.df) {
  names(.df) <- gsub("Back-.*", "estimate", names(.df))
  names(.df) <- gsub("Estimate", "model.est", names(.df))
  names(.df) <- gsub(".*RSE.*", "rse", names(.df))
  names(.df) <- gsub("SE", "std.error", names(.df))
  names(.df) <- gsub("CI Lower", "conf.low", names(.df))
  names(.df) <- gsub("CI Upper", "conf.high", names(.df))
  names(.df) <- gsub("BSV.*", "bsv", names(.df))
  names(.df) <- gsub("Shrink.*", "shrink", names(.df))
  return(.df)
}

##' @export
confint.nlmixr2FitCore <- function(object, parm, level = 0.95, ...) {
  .extra <- list(...)
  .exponentiate <- ifelse(any(names(.extra) == "exponentiate"), .extra$exponentiate, FALSE)
  .ciNames <- ifelse(any(names(.extra) == "ciNames"), .extra$ciNames, TRUE)
  .df <- .fixNames(object$parFixedDf)
  if (!missing(parm)) .df <- .df[parm, ]
  .zv <- qnorm(1 - (1 - level) / 2)
  .low <- .df$model.est - .df$std.error * .zv
  .hi <- .df$model.est + .df$std.error * .zv
  if (is.na(.exponentiate)) {
    .exp <- abs(exp(.df$model.est) - .df$estimate) < 1e-6
    .hi[.exp] <- exp(.hi[.exp])
    .low[.exp] <- exp(.low[.exp])
  } else if (.exponentiate) {
    .hi <- exp(.hi)
    .low <- exp(.low)
  }
  .df <- data.frame(model.est = .df$model.est, estimate = .df$estimate, conf.low = .low, conf.high = .hi, row.names=rownames(.df))
  if (.ciNames) names(.df)[3:4] <- paste(c((1 - level) / 2, (1 - (1 - level) / 2)) * 100, "%")
  .df
}

##' @export
confint.nlmixr2FitCoreSilent <- confint.nlmixr2FitCore

.nlmixr2TidyFixed <- function(x, ..., .ranpar = FALSE) {
  rxode2::rxReq("tibble")
  .extra <- list(...)
  .conf.int <- ifelse(any(names(.extra) == "conf.int"), .extra$conf.int, ifelse(any(names(.extra) == "conf.level"), TRUE, FALSE))
  .conf.level <- ifelse(any(names(.extra) == "conf.level"), .extra$conf.level, 0.95)
  .exponentiate <- ifelse(any(names(.extra) == "exponentiate"), .extra$exponentiate, FALSE)
  .quick <- ifelse(any(names(.extra) == "quick"), .extra$quick, FALSE)
  .rse <- ifelse(any(names(.extra) == "rse"), .extra$rse, FALSE)
  .bsv <- ifelse(any(names(.extra) == "bsv"), .extra$bsv, FALSE)
  .shrink <- ifelse(any(names(.extra) == "shrink"), .extra$shrink, FALSE)
  if (.quick) warning("quick does not do anything for nlmixr2 fit objects")
  .df <- .fixNames(x$parFixedDf)
  .exp <- abs(exp(.df$model.est) - .df$estimate) < 1e-6
  if (is.na(.exponentiate)) {
    ## se(exp(x)) ~= exp(mu_x)*se_x
    ##
    ## Norris N. The Standard Errors of the Geometric and Harmonic Means and Their Application to Index Numbers
    ## Ann. Math. Statist. Volume 11, Number 4 (1940), 445-448.
    .df$std.error[.exp] <- exp(.df$model.est[.exp]) * .df$std.error[.exp]
  } else if (.exponentiate) {
    .df$std.error <- exp(.df$model.est) * .df$std.error
    .df$estimate <- exp(.df$model.est)
  } else if (!.exponentiate) {
    .df$estimate <- .df$model.est
  }
  .df$statistic <- .df$estimate / .df$std.error
  .df$p.value <- stats::pt(.df$statistic, nobs(x) - attr(logLik(x), "df"), lower.tail = FALSE)
  if (.conf.int) {
    .ci <- confint.nlmixr2FitCore(x, level = .conf.level, ciNames = FALSE, exponentiate = .exponentiate)
    .df$conf.low <- .ci$conf.low
    .df$conf.high <- .ci$conf.high
  } else {
    .df <- .df[, regexpr("^conf[.]", names(.df)) == -1]
  }
  if (!.rse) {
    .df <- .df[, names(.df) != "rse"]
  }
  if (!.bsv) {
    .df <- .df[, names(.df) != "bsv"]
  }
  if (!.shrink) {
    .df <- .df[, names(.df) != "shrink"]
  }
  .ini <- x$ui$iniDf$err[!is.na(x$ui$iniDf$ntheta)]
  ## effect   group   term            estimate std.error statistic
  .df <- data.frame(
    effect = "fixed",
    term = row.names(.df), .df, stringsAsFactors = FALSE
  )
  if (!.ranpar) {
    .df <- .df[is.na(x$ui$iniDf$err[!is.na(x$ui$iniDf$ntheta)]), ]
  } else {
    .tmp <- x$ui$iniDf$err[!is.na(x$ui$iniDf$ntheta)]
    .df <- .df[!is.na(.tmp), ]
    .tmp <- .tmp[!is.na(.tmp)]
    .df$group <- paste0("Residual(", .tmp, ")")
    .df$effect <- "ran_pars"
  }
  tibble::as_tibble(.df)
}

.nlmixr2TidyRandom <- function(x, ...) {
  rxode2::rxReq("tibble")
  .d <- dim(x$omegaR)
  if (.d[1] > 0) {
    .tmp <- stack(x$eta[, -1])
    .df <- data.frame(group = "ID", level = x$eta$ID, term = .tmp$ind, estimate = .tmp$values)
    return(tibble::as_tibble(.df))
  } else {
    return(NULL)
  }
}

.nlmixr2TidyRandomPar <- function(x, ...) {
  rxode2::rxReq("tibble")
  .pars <- .getR(x$omegaR, TRUE)
  if (length(.pars) > 0) {
    .p1 <- data.frame(
      effect = "ran_pars", group = "ID", term = names(.pars), estimate = .pars, std.error = NA_real_,
      statistic = NA_real_, p.value = NA_real_, stringsAsFactors = FALSE
    ) %>%
      .reorderCols()
    .p2 <- data.frame(.nlmixr2TidyFixed(x, .ranpar = TRUE), stringsAsFactors = FALSE) %>%
      .reorderCols()
    .df <- rbind(.p1, .p2)
    for (.v in c("statistic", "p.value", "std.error")) {
      if (all(is.na(.df[, .v]))) {
        .df <- .df[, names(.df) != .v]
      }
    }
    return(tibble::as_tibble(.df))
  } else {
    return(NULL)
  }
}
##   effect   group    term                  estimate std.error statistic
##   <chr>    <chr>    <chr>                    <dbl>     <dbl>     <dbl>
## 1 fixed    NA       (Intercept)           251.          6.82     36.8
## 2 fixed    NA       Days                   10.5         1.55      6.77
## 3 ran_pars Subject  sd__(Intercept)        24.7        NA        NA
## 4 ran_pars Subject  sd__Days                5.92       NA        NA
## 5 ran_pars Subject  cor__(Intercept).Days   0.0656     NA        NA
## 6 ran_pars Residual sd__Observation        25.6        NA        NA

## Row names and order taken & adapted from
## https://github.com/bbolker/broom.mixed/blob/master/R/utilities.R#L238-L248
.reorderCols <- function(x) {
  allCols <- c(
    "response", "effect",
    "component", ## glmmTMB, brms
    "group", "level", "term", "index", "estimate",
    "std.error", "statistic",
    "df", "p.value",
    "conf.low", "conf.high", "rhat", "ess"
  )
  return(x[, intersect(allCols, names(x))])
}

.coefPar <- function(x, exponentiate = FALSE, ...) {
  rxode2::rxReq("tibble")
  .d <- dim(x$omegaR)
  if (.d[1] == 0) {
    return(NULL)
  }
  .muRef <- x$ui$muRefDataFrame
  .curEval <- x$ui$muRefCurEval
  .theta <- x$theta
  .df <- .fixNames(x$parFixedDf)
  if (is.na(exponentiate)) {
    .exp <- abs(exp(.df$model.est) - .df$estimate) < 1e-6
  } else if (exponentiate) {
    .exp <- setNames(rep(TRUE, length(.theta)), names(.theta))
  } else {
    .exp <- setNames(rep(FALSE, length(.theta)), names(.theta))
  }
  .eta <- x$eta
  .noMuRef <- c()
  .x <- setNames(
    data.frame(lapply(names(.eta), function(eta) {
      .w <- which(.muRef$eta == eta)
      if (length(.w) == 1L) {
        .thetaName <- .muRef$theta[.w]
        .ret <- .eta[[eta]] + .theta[.thetaName]
        if (any(.thetaName == names(.exp))) {
          if (.exp[.thetaName]) {
            .ret <- exp(.ret)
          }
        }
      } else {
        if (eta != "ID") {
          warning(sprintf("the parameter '%s' is not mu-referenced and the coef will not be returned", eta), call.=FALSE)
          .noMuRef <<- c(.noMuRef, eta)
        }
        .ret <- .eta[[eta]]
      }
      return(.ret)
    })), sapply(names(.eta), function(eta) {
      .w <- which(.muRef$eta == eta)
      if (length(.w) == 1L) {
        .thetaName <- .muRef$theta[.w]
        return(.thetaName)
      }
      return(eta)
    })
  )
  .tmp <- stack(.x[, -1])
  .df <- data.frame(group = "ID", level = .x$ID, term = .tmp$ind, estimate = .tmp$values)
  return(tibble::as_tibble(.df))
}

tidy.nlmixr2FitCore <- function(x, ...) {
  rxode2::rxReq("tibble")
  rxode2::rxReq("dplyr")
  .extra <- list(...)
  if (any(names(.extra) == "effects")) {
    .effects <- .extra$effects
  } else {
    if (any(names(.extra) == "effect")) {
      .effects <- .extra$effect
    } else {
      .effects <- c("fixed", "ran_pars")
    }
  }
  .effects <- match.arg(.effects, c("fixed", "random", "ran_vals", "ran_pars", "ran_coef"),
    several.ok = TRUE
  )
  .ret <- list()
  if (any(.effects == "fixed")) {
    .ret$fixed <- .nlmixr2TidyFixed(x, ...)
  }
  if (any(.effects == "random") || any(.effects == "ran_vals")) {
    .ret$ran_vals <- .nlmixr2TidyRandom(x, ...)
  }
  if (any(.effects == "ran_pars")) {
    .ret$ran_pars <- .nlmixr2TidyRandomPar(x, ...)
  }
  if (any(.effects == "ran_coef")) {
    .ret$ran_coef <- .coefPar(x, ...)
  }
  if (all(unlist(lapply(seq_along(.ret), is.null)))) {
    return(NULL)
  }
  return(dplyr::bind_rows(.ret, .id = "effect") %>%
    tibble::as_tibble() %>%
    .reorderCols())
}

tidy.nlmixr2FitCoreSilent <- tidy.nlmixr2FitCore

glance.nlmixr2FitCore <- function(x, ...) {
  rxode2::rxReq("tibble")
  .lst <- list(...)
  if (any(names(.lst) == "type")) {
    setOfv(x, type = .lst$type)
  }
  .aic <- AIC(x) ## To calculate AIC if needed
  .df <- x$objDf
  if (length(x$logLik) > 1L) {
    .df <- .df[which(tolower(rownames(x$objDf)) == tolower(x$ofvType)), ]
  }
  names(.df) <- gsub("Log-likelihood", "logLik", names(.df))
  names(.df) <- gsub("Condition#(Cor)", "conditionNumberCor", names(.df))
  names(.df) <- gsub("Condition#(Cov)", "conditionNumberCov", names(.df))
  tibble::as_tibble(.df)
}

glance.nlmixr2FitCoreSilent <- glance.nlmixr2FitCore

augment.nlmixr2FitCore <- function(x, ...) {
  stop("augment is not yet implemented for nlmixr2 models",
       call.=FALSE)
}

augment.nlmixr2FitCoreSilent <- augment.nlmixr2FitCore
