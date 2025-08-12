#' Control options foir the FOI estimation method
#'
#' This is related to the focei methods and uses most of their control
#' options.  Some are ignored, `posthoc` is an extra parameter
#'
#' @inheritParams foceiControl
#' @param ... Parameters used in the default `foceiConrol()`
#' @param interaction Interaction term for the model; ignored by fo
#' @param fo Logical indicating if the estimation method is FO (first
#'   order), but this is controlled by the estimation method so this
#'   is ignored.
#' @param posthoc Logical indicating if the estimation method should
#'   calculate `foce` posthoc predicted parameters.
#' @return foiControl object
#' @export
#' @author Matthew L. Fidler
#' @examples
#'
#' foiControl()
foiControl <- function(sigdig=3,
                       ...,
                       posthoc=TRUE,
                       interaction=NULL,
                       fo=NULL) {
  checkmate::assertLogical(posthoc, len=1, any.missing=FALSE, null.ok=FALSE)
  .control <- foceiControl(sigdig=sigdig, ..., interaction=1L, fo=TRUE)
  class(.control) <- NULL
  .control <- c(.control, list(posthoc=posthoc))
  class(.control) <- "foiControl"
  .control
}

#' @rdname nmObjHandleControlObject
#' @export
nmObjHandleControlObject.foiControl <- function(control, env) {
  assign("foiControl", control, envir=env)
}

#' @rdname getValidNlmixrControl
#' @export
getValidNlmixrCtl.foi <- function(control) {
  .ctl <- control[[1]]
  .cls <- class(control)[1]
  if (is.null(.ctl)) .ctl <- foiControl()
  if (is.null(attr(.ctl, "class")) && is(.ctl, "list"))
    .ctl <- do.call("foiControl", .ctl)
  if (inherits(.ctl, "foceiControl") ||
        inherits(.ctl, "foceControl") ||
        inherits(.ctl, "foControl")) {
    .minfo(paste0("converting ", class(.ctl)[1], " to foiControl"))
    class(.ctl) <- NULL
    .ctl <- do.call(foiControl, .ctl)
  } else if (inherits(.ctl, "foceControl")) {
    .minfo(paste0("converting foceControl to foiControl"))
    class(.ctl) <- NULL
    .ctl <- do.call(foiControl, .ctl)
  } else if (!inherits(.ctl, "foiControl")) {
    .minfo(paste0("invalid control for `est=\"", .cls, "\"`, using default"))
    .ctl <- foiControl()
  } else {
    .ctl <- do.call(foiControl, .ctl)
  }
  .ctl
}


#' @rdname nmObjGetControl
#' @export
nmObjGetControl.foi <- function(x, ...) {
  .env <- x[[1]]
  for (.name in c("foiControl", "control",
                  "foControl", "foceControl",
                  "foceiControl", "foceiControl0")) {
    if (exists(.name, .env)) {
      .control <- get(.name, .env)
      if (inherits(.control, "foiControl")) return(.control)
      .ret <- try(suppressMessages(getValidNlmixrCtl.foi(list(.control))), silent=TRUE)
      if (inherits(.ret, "foiControl")) return(.ret)
    }
  }
  stop("cannot find foi related control object", call.=FALSE)
}


.foiControlToFoceiControl <- function(env, assign=TRUE) {
  .foiControl <- env$foiControl
  .ui <- env$ui
  .n <- names(.foiControl)
  .w <- which(.n == "posthoc")
  .n <- .n[-.w]
  .foceiControl <- setNames(lapply(.n,
                                   function(n) {
                                     if (n == "maxInnerIterations" &&
                                           !.foiControl$posthoc) {
                                       return(0L)
                                     }
                                     .foiControl[[n]]
                                   }), .n)
  class(.foceiControl) <- "foceiControl"
  if (assign) env$control <- .foceiControl
  .foceiControl
}


#'@rdname nlmixr2Est
#'@export
nlmixr2Est.foi <- function(env, ...) {
  .ui <- env$ui
  rxode2::assertRxUiTransformNormal(.ui, " for the estimation routine 'foi'", .var.name=.ui$modelName)
  rxode2::assertRxUiRandomOnIdOnly(.ui, " for the estimation routine 'foi'", .var.name=.ui$modelName)
  rxode2::assertRxUiMixedOnly(.ui, " for the estimation routine 'foi'", .var.name=.ui$modelName)
  .foceiFamilyControl(env, ..., type="foiControl")
  .control <- .ui$control
  .posthoc <- .control$posthoc
  rxode2::rxAssignControlValue(.ui, "interaction", 0L)
  rxode2::rxAssignControlValue(.ui, "covMethod", 0L)
  rxode2::rxAssignControlValue(.ui, "fo", TRUE)
  rxode2::rxAssignControlValue(.ui, "boundTol", 0)
  rxode2::rxAssignControlValue(.ui, "compress", 0L)
  on.exit({
    if (exists("control", envir=.ui)) {
      rm("control", envir=.ui)
    }
  })
  env$skipTable <- TRUE
  .ret <- .foceiFamilyReturn(env, .ui, ...)
  .objDf <- .ret$objDf
  .ui <- .ret$ui # This has the final estimates

  ## Now the posthoc/table step
  env$foiControl <- .control
  .foceiFamilyControl(env, ..., type="foiControl")
  .foiControlToFoceiControl(env)
  .ui$control <- env$control
  rxode2::rxAssignControlValue(.ui, "interaction", 1L)
  rxode2::rxAssignControlValue(.ui, "covMethod", .control$covMethod)
  rxode2::rxAssignControlValue(.ui, "fo", FALSE)
  rxode2::rxAssignControlValue(.ui, "boundTol", .control$boundTol)
  rxode2::rxAssignControlValue(.ui, "compress", .control$compress)
  if (.control$posthoc) {
    rxode2::rxAssignControlValue(.ui, "maxInnerIterations", .control$maxInnerIterations)
  } else {
    rxode2::rxAssignControlValue(.ui, "maxInnerIterations", 0L)
  }
  rm(list="skipTable", envir=env)

  rxode2::rxAssignControlValue(.ui, "maxOuterIterations", 0L)
  rxode2::rxAssignControlValue(.ui, "calcTables", .control$calcTables)
  env$control <- .control
  .ret <- .foceiFamilyReturn(env, .ui, ..., method="FOI", est="foi")
  .ret$est <- "foi"
  assign("foiControl", .control, envir=.ret$env)
  assign("control", env$control, envir=.ret$env)
  rm("control", envir=.ret$env)
  .addObjDfToReturn(.ret, .objDf)
  .ret
}
attr(nlmixr2Est.foi, "covPresent") <- TRUE
