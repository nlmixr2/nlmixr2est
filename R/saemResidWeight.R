#' Functions a known residual weight may call
#'
#' Plain R evaluates them the same way rxode2 does, so the weight can be
#' computed from the data in R.
#' @noRd
.saemResidWeightFuns <- c(
  "+",
  "-",
  "*",
  "/",
  "^",
  "(",
  "exp",
  "log",
  "sqrt",
  "abs",
  "log10",
  "log2",
  "log1p",
  "expm1"
)

#' Split an error argument into a theta times a known weight
#'
#' @param expr right hand side defining the error argument
#' @param theta the theta to factor out
#' @return the weight expression, or `NULL` when `expr` is not `theta * w`
#' @noRd
.saemFactorTheta <- function(expr, theta) {
  if (is.name(expr)) {
    return(if (identical(as.character(expr), theta)) 1 else NULL)
  }
  if (!is.call(expr)) {
    return(NULL)
  }
  .op <- as.character(expr[[1]])
  if (.op == "(" && length(expr) == 2L) {
    return(.saemFactorTheta(expr[[2]], theta))
  }
  if (!(.op %in% c("*", "/")) || length(expr) != 3L) {
    return(NULL)
  }
  .inL <- theta %in% all.vars(expr[[2]])
  .inR <- theta %in% all.vars(expr[[3]])
  if (.inL == .inR) {
    return(NULL)
  }
  if (.inL) {
    .w <- .saemFactorTheta(expr[[2]], theta)
    if (is.null(.w)) {
      return(NULL)
    }
    if (identical(.w, 1) && .op == "*") {
      return(expr[[3]])
    }
    return(call(.op, .w, expr[[3]]))
  }
  # theta in a divisor is not a scale factor
  if (.op == "/") {
    return(NULL)
  }
  .w <- .saemFactorTheta(expr[[3]], theta)
  if (is.null(.w)) {
    return(NULL)
  }
  if (identical(.w, 1)) {
    return(expr[[2]])
  }
  call("*", expr[[2]], .w)
}

#' Is this call only built from known functions?
#'
#' @param expr expression
#' @return logical
#' @noRd
.saemResidWeightFunsOk <- function(expr) {
  if (!is.call(expr)) {
    return(TRUE)
  }
  if (!is.name(expr[[1]]) || !(as.character(expr[[1]]) %in% .saemResidWeightFuns)) {
    return(FALSE)
  }
  all(vapply(as.list(expr[-1]), .saemResidWeightFunsOk, logical(1)))
}

#' Is a model line an assignment to `name`?
#'
#' @param e model line
#' @param name variable name
#' @return logical
#' @noRd
.saemIsAssignTo <- function(e, name) {
  is.call(e) &&
    length(e) == 3L &&
    (identical(e[[1]], as.name("<-")) || identical(e[[1]], as.name("=")) || identical(e[[1]], as.name("~"))) &&
    identical(e[[2]], as.name(name))
}

#' Additive or proportional residuals that are a theta times a known weight
#'
#' Finds normal endpoints written as `add(theta * w)` or `prop(theta * w)`,
#' where `w` only depends on `time` and covariates.  saem fits these on its
#' closed-form residual path: an additive one by scaling the observations and
#' predictions by `1/w`, a proportional one with `w` as a per-observation
#' multiplier of the residual SD.
#'
#' @param ui rxode2 ui
#' @param conds endpoint conditions with a modeled residual
#'   (`.saemModeledResidualCond()`)
#' @return list of specs (`cond`, `cmt`, `type`, `theta`, `w`, `def`, `err`),
#'   one per qualifying endpoint
#' @noRd
.saemResidWeightSpecs <- function(ui, conds) {
  .pred <- ui$predDf
  .iniDf <- ui$iniDf
  .lst <- ui$lstExpr
  .thetas <- .iniDf$name[!is.na(.iniDf$ntheta) & is.na(.iniDf$err)]
  .known <- c("time", "t", ui$allCovs)
  .ret <- list()
  for (.cond in conds) {
    .i <- which(.pred$cond == .cond)
    if (length(.i) != 1L) {
      next
    }
    .p <- .pred[.i, , drop = FALSE]
    .type <- paste(.p$errType)
    if (
      !identical(paste(.p$distribution), "norm") ||
        !(.type %in% c("add", "prop")) ||
        !identical(paste(.p$transform), "untransformed")
    ) {
      next
    }
    .col <- if (.type == "add") "a" else "b"
    .others <- setdiff(intersect(c("a", "b", "c", "d", "e", "f", "lambda"), names(.pred)), .col)
    if (any(!is.na(unlist(.p[, .others, drop = TRUE])))) {
      next
    }
    .a <- paste(.p[[.col]])
    .errLine <- .p$line
    .def <- which(vapply(.lst, .saemIsAssignTo, logical(1), name = .a))
    if (length(.def) != 1L || .def >= .errLine) {
      next
    }
    .rhs <- .lst[[.def]][[3]]
    .th <- intersect(all.vars(.rhs), .thetas)
    if (length(.th) != 1L) {
      next
    }
    # the residual theta and argument may not be used anywhere else
    .rest <- .lst[-c(.def, .errLine)]
    .used <- unique(unlist(lapply(.rest, all.vars)))
    if (any(c(.th, .a) %in% .used)) {
      next
    }
    .w <- .saemFactorTheta(.rhs, .th)
    if (is.null(.w) || !all(all.vars(.w) %in% .known) || !.saemResidWeightFunsOk(.w)) {
      next
    }
    .ret[[length(.ret) + 1L]] <- list(
      cond = .cond,
      cmt = .p$cmt,
      type = .type,
      theta = .th,
      w = .w,
      def = .lst[[.def]],
      err = .lst[[.errLine]]
    )
  }
  .ret
}

#' Name of the line that holds a weighted residual's place in the fit's model
#'
#' The scaled prediction saem fits for an additive residual, the (unused)
#' weight for a proportional one.
#'
#' @param spec one `.saemResidWeightSpecs()` element
#' @return variable name
#' @noRd
.saemResidWeightVar <- function(spec) {
  paste0("rx.", spec$cond, if (identical(spec$type, "prop")) ".w" else ".wpred")
}

#' Rebuild a ui from changed model lines, keeping its ini block
#'
#' @param ui rxode2 ui
#' @param lst new model lines
#' @return rxode2 ui (decompressed)
#' @noRd
.saemResidWeightRebuild <- function(ui, lst) {
  ui <- rxode2::rxUiDecompress(ui)
  .model <- as.call(c(list(as.name("model")), list(as.call(c(list(as.name("{")), lst)))))
  .ini <- as.expression(lotri::as.lotri(ui$iniDf))
  .ini[[1]] <- quote(`ini`)
  .fun <- .getUiFunFromIniAndModel(ui, .ini, .model)
  .newUi <- rxode2::rxUiDecompress(.fun())
  assign("modelName", ui$modelName, envir = .newUi)
  for (.v in c("boundedTransforms", "control")) {
    if (exists(.v, envir = ui, inherits = FALSE)) {
      assign(.v, get(.v, envir = ui), envir = .newUi)
    }
  }
  .newUi
}

#' Fit a weighted residual as a plain additive or proportional one
#'
#' `cp ~ add(add.sd * w)` becomes `rx.cp.wpred <- cp / (w)` and
#' `rx.cp.wpred ~ add(add.sd) | cp`; saem divides the endpoint's `DV` and
#' `LIMIT` by the same `w` (`.saemResidWeightData()`), which leaves the
#' likelihood unchanged up to a constant.  `cp ~ prop(prop.sd * w)` becomes
#' `cp ~ prop(prop.sd)`, with `w` passed to saem as a per-observation
#' multiplier of the residual SD.
#'
#' @param ui rxode2 ui
#' @param specs `.saemResidWeightSpecs()`
#' @return rewritten ui
#' @noRd
.saemResidWeightUi <- function(ui, specs) {
  .lst <- ui$lstExpr
  for (.s in specs) {
    .defI <- which(vapply(.lst, identical, logical(1), y = .s$def))
    .errI <- which(vapply(.lst, identical, logical(1), y = .s$err))
    .var <- as.name(.saemResidWeightVar(.s))
    .theta <- as.name(.s$theta)
    if (identical(.s$type, "prop")) {
      .lst[[.defI]] <- bquote(.(.var) <- .(.s$w))
      .lst[[.errI]] <- bquote(.(as.name(.s$cond)) ~ prop(.(.theta)))
    } else {
      .lst[[.defI]] <- bquote(.(.var) <- .(as.name(.s$cond)) / (.(.s$w)))
      .lst[[.errI]] <- bquote(.(.var) ~ add(.(.theta)) | .(as.name(.s$cond)))
    }
  }
  .saemResidWeightRebuild(ui, .lst)
}

#' Put a weighted residual back the way the model wrote it
#'
#' @param ui the fit's ui (`.saemResidWeightUi()` form)
#' @param specs `.saemResidWeightSpecs()`
#' @return ui with the original residual lines
#' @noRd
.saemResidWeightRestoreUi <- function(ui, specs) {
  if (length(specs) == 0L) {
    return(ui)
  }
  .lst <- ui$lstExpr
  for (.s in specs) {
    .var <- .saemResidWeightVar(.s)
    .errVar <- if (identical(.s$type, "prop")) .s$cond else .var
    .defI <- which(vapply(
      .lst,
      function(e) {
        .saemIsAssignTo(e, .var) && identical(e[[1]], as.name("<-"))
      },
      logical(1)
    ))
    .errI <- which(vapply(
      .lst,
      function(e) {
        is.call(e) && identical(e[[1]], as.name("~")) && identical(e[[2]], as.name(.errVar))
      },
      logical(1)
    ))
    if (length(.defI) != 1L || length(.errI) != 1L) {
      next
    }
    .lst[[.defI]] <- .s$def
    .lst[[.errI]] <- .s$err
  }
  .ui <- .saemResidWeightRebuild(ui, .lst)
  .ui
}

#' Apply the known residual weights to saem's data
#'
#' An additive endpoint's observations are divided by the weight; a
#' proportional endpoint's weight becomes a per-observation residual SD
#' multiplier.
#'
#' @param data saem's processed data (`$dataSav`)
#' @param specs `.saemResidWeightSpecs()`
#' @param nendpnt number of endpoints in the model
#' @return list with the scaled `data`, `logW` (the sum of `log(w)` over the
#'   scaled observations) and `resW` (the residual SD multiplier of each
#'   observation, `NULL` without a proportional weight)
#' @noRd
.saemResidWeightData <- function(data, specs, nendpnt) {
  .logW <- 0
  .resW <- NULL
  if (length(specs) == 0L) {
    return(list(data = data, logW = .logW, resW = .resW))
  }
  .names <- names(data)
  .time <- data[[which(toupper(.names) == "TIME")[1]]]
  .evid <- data[[which(toupper(.names) == "EVID")[1]]]
  .cmtCol <- which(toupper(.names) == "CMT")
  .dv <- which(toupper(.names) == "DV")[1]
  .lim <- which(toupper(.names) == "LIMIT")
  .env <- c(as.list(data), list(time = .time, t = .time))
  for (.s in specs) {
    .rows <- .evid == 0
    if (nendpnt > 1L && length(.cmtCol) == 1L) {
      .rows <- .rows & data[[.cmtCol]] == .s$cmt
    }
    .w <- eval(.s$w, envir = .env, enclos = baseenv())
    .w <- rep_len(as.numeric(.w), nrow(data))[.rows]
    if (any(!is.finite(.w) | .w <= 0)) {
      stop(sprintf("residual weight for '%s' must be positive and finite", .s$cond), call. = FALSE)
    }
    if (identical(.s$type, "prop")) {
      if (is.null(.resW)) {
        .resW <- rep(1, sum(.evid == 0))
      }
      .resW[.rows[.evid == 0]] <- .w
      next
    }
    data[[.dv]][.rows] <- data[[.dv]][.rows] / .w
    if (length(.lim) == 1L) {
      .l <- data[[.lim]][.rows]
      .fin <- is.finite(.l)
      .l[.fin] <- .l[.fin] / .w[.fin]
      data[[.lim]][.rows] <- .l
    }
    .logW <- .logW + sum(log(.w))
  }
  list(data = data, logW = .logW, resW = .resW)
}
