#' Handle Model Object
#'
#' @param model model list should have at least:
#'
#' - `predOnly` -- this is the prediction model with all the left
#'    handed equations added so they will be added the table.  The
#'    model should have `rx_pred_`, the model based prediction, as the
#'    first defined lhs component.  The second component should be
#'    `rx_r_`, the variance of the prediction.  These variables may
#'    change based on distribution type.  In additional all
#'    interesting calculated variables should be included.
#'
#' - `predNoLhs` -- This is the prediction model.  It only has the
#'    prediction and no left handed equations.
#'
#' @return This returns the `$model` object for a fit.  It is a s3
#'   method because it may be different between different model types
#'
#' @param env Environment for the fit information
nmObjHandleModelObject <- function(model, env) {
  on.exit({
    if (exists("model", envir=env)){
      rm("model", envir=env)
    }})
  UseMethod("nmObjHandleModelObject")
}

#' @rdname nmObjHandleModelObject
#' @export
nmObjHandleModelObject.saemModelList <- function(model, env) {
  assign("saemModel", model, envir=env)
}

 #' @rdname nmObjHandleModelObject
#' @export
nmObjHandleModelObject.foceiModelList <- function(model, env) {
  assign("foceiModel", model, envir=env)
}

#' @rdname nmObjHandleModelObject
#' @export
nmObjHandleModelObject.default <- function(model, env) {
  stop("cannot figure out how to handle the model, add method for `nmObjHandleModelObject.", class(model)[1], "`",
       call.=FALSE)
}

#' Handle the control object
#'
#'
#' @param control Control object
#' @param env fit environment
#' @return Nothing, called for side effects
#' @author Matthew L. Fidler
#' @export
nmObjHandleControlObject <- function(control, env) {
  on.exit({
    if (exists("control", envir=env)) {
      rm("control", envir=env)
    }
  })
  UseMethod("nmObjHandleControlObject")
}

#' @rdname nmObjHandleControlObject
#' @export
nmObjHandleControlObject.foceiControl <- function(control, env) {
  assign("foceiControl0", control, envir=env)
}


#' @rdname nmObjHandleControlObject
#' @export
nmObjHandleControlObject.saemControl <- function(control, env) {

  assign("saemControl", control, envir=env)
}


#' @rdname nmObjHandleControlObject
#' @export
nmObjHandleControlObject.default <- function(control, env) {
  stop("cannot figure out how to handle the model, add method for `nmObjHandleControlObject.", class(control)[1], "`",
       call.=FALSE)
}

#' Method for getting focei compatible control object from nlmixr object
#'
#' @param x nlmixr composed fit object
#'
#' @param ... Other parameters
#'
#' @return foceiControl translated from current control
#'
#' @export
nmObjGetFoceiControl <- function(x, ...) {
  UseMethod("nmObjGetFoceiControl")
}

#' @rdname nmObjGetFoceiControl
#' @export
nmObjGetFoceiControl.default <- function(x, ...) {
  .env <- x[[1]]
  if (exists("foceiControl0", .env)) {
    return(get("foceiControl0", .env))
  } else {
    stop("cannot figure out how to make/retrieve the focei control\nmissing 'nmObjGetFoceiControl.",
         class(x)[1], "'",
         call.=FALSE)
  }
}

#' Get control object from fit
#'
#' @param x nlmixr fit object
#' @param ... Other parameters
#' @return Control object of estimation method
#' @author Matthew L. Fidler
#' @export
nmObjGetControl <- function(x, ...) {
  UseMethod("nmObjGetControl")
}

#' @rdname nmObjGetControl
#' @export
nmObjGetControl.focei <- function(x, ...) {
  .env <- x[[1]]
  if (exists("foceiControl0", .env)) {
    .control <- get("foceiControl0", .env)
    if (inherits(.control, "foceiControl")) return(.control)
  }
  if (exists("control", .env)) {
    .control <- get("control", .env)
    if (inherits(.control, "foceiControl")) return(.control)
  }
  stop("cannot find focei related control object", call.=FALSE)
}


#' @rdname nmObjGetControl
#' @export
nmObjGetControl.saem <- function(x, ...) {
  .env <- x[[1]]
  if (exists("saemControl", .env)) {
    .control <- get("saemControl", .env)
    if (inherits(.control, "saemControl")) return(.control)
  }
  if (exists("control", .env)) {
    .control <- get("control", .env)
    if (inherits(.control, "saemControl")) return(.control)
  }
  stop("cannot find saem related control object", call.=FALSE)
}

#' @rdname nmObjGetControl
#' @export
nmObjGetControl.default <- function(x, ...) {
  ## stop("cannot figure out get the control, add method for `nmObjHandleControlObject.", class(x)[1], "`",
  ##      call.=FALSE)
  NULL
}
