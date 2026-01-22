#' Control options for the FO estimation method
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
#' @return foControl object
#' @export
#' @author Matthew L. Fidler
#' @examples
#'
#' foControl()
foControl <- function(sigdig=3,
                      ...,
                      posthoc=TRUE,
                      interaction=NULL,
                      fo=NULL) {
  checkmate::assertLogical(posthoc, len=1, any.missing=FALSE, null.ok=FALSE)
  .control <- foceiControl(sigdig=sigdig, ..., interaction=0L, fo=TRUE)
  class(.control) <- NULL
  .control <- c(.control, list(posthoc=posthoc))
  class(.control) <- "foControl"
  .control
}

#' @rdname nmObjHandleControlObject
#' @export
nmObjHandleControlObject.foControl <- function(control, env) {
  ## eval(rxode2::rxUiDeparse(control, "control"))
  assign("foControl", control, envir=env)
}

#' @rdname getValidNlmixrControl
#' @export
getValidNlmixrCtl.fo <- function(control) {
  .ctl <- control[[1]]
  .cls <- class(control)[1]
  if (is.null(.ctl)) .ctl <- foControl()
  if (is.null(attr(.ctl, "class")) && is(.ctl, "list"))
    .ctl <- do.call("foControl", .ctl)
  if (inherits(.ctl, "foceiControl") ||
        inherits(.ctl, "foceControl") ||
        inherits(.ctl, "foiControl")) {
    .minfo(paste0("converting ", class(.ctl)[1], " to foControl"))
    class(.ctl) <- NULL
    .ctl <- do.call(foControl, .ctl)
  } else if (inherits(.ctl, "foceControl")) {
    .minfo(paste0("converting foceControl to foControl"))
    class(.ctl) <- NULL
    .ctl <- do.call(foControl, .ctl)
  } else if (!inherits(.ctl, "foControl")) {
    .minfo(paste0("invalid control for `est=\"", .cls, "\"`, using default"))
    .ctl <- foControl()
  } else {
    .ctl <- do.call(foControl, .ctl)
  }
  .ctl
}


#' @rdname nmObjGetControl
#' @export
nmObjGetControl.fo <- function(x, ...) {
  .env <- x[[1]]
  for (.name in c("foControl", "control",
                  "foiControl", "foceControl",
                  "foceiControl", "foceiControl0")) {
    if (exists(.name, .env)) {
      .control <- get(.name, .env)
      if (inherits(.control, "foControl")) return(.control)
      .ret <- try(suppressMessages(getValidNlmixrCtl.fo(list(.control))), silent=TRUE)
      if (inherits(.ret, "foControl")) return(.ret)
    }
  }
  stop("cannot find fo related control object", call.=FALSE)
}


.foControlToFoceiControl <- function(env, assign=TRUE) {
  .foControl <- env$foControl
  .ui <- env$ui
  .n <- names(.foControl)
  .w <- which(.n == "posthoc")
  .n <- .n[-.w]
  .foceiControl <- setNames(lapply(.n,
         function(n) {
           if (n == "maxInnerIterations" && !.foControl$posthoc) {
             return(0L)
           }
          .foControl[[n]]
         }), .n)
  class(.foceiControl) <- "foceiControl"
  if (assign) env$control <- .foceiControl
  .foceiControl
}


#'@rdname nlmixr2Est
#'@export
nlmixr2Est.fo <- function(env, ...) {
  .ui <- env$ui
  rxode2::assertRxUiTransformNormal(.ui, " for the estimation routine 'fo'", .var.name=.ui$modelName)
  rxode2::assertRxUiRandomOnIdOnly(.ui, " for the estimation routine 'fo'", .var.name=.ui$modelName)
  rxode2::assertRxUiMixedOnly(.ui, " for the estimation routine 'fo'", .var.name=.ui$modelName)
  .foceiFamilyControl(env, ..., type="foControl")
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
  env$foControl <- .control
  .foceiFamilyControl(env, ..., type="foControl")
  .foControlToFoceiControl(env)
  .ui$control <- env$control
  rxode2::rxAssignControlValue(.ui, "interaction", 0L)
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
  .ret <- .foceiFamilyReturn(env, .ui, ..., method="FO", est="fo")
  .ret$est <- "fo"
  assign("foControl", .control, envir=.ret$env)
  assign("control", env$control, envir=.ret$env)
  rm("control", envir=.ret$env)
  .addObjDfToReturn(.ret, .objDf)
  .ret
}
attr(nlmixr2Est.fo, "covPresent") <- TRUE
