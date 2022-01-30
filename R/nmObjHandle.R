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
#' @param env Environment for the fit information
nmObjHandleModelObject <- function(model, env) {
  on.exit(rm("model", envir=env))
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
nmObjHandleControlObject <- function(control, env) {
  on.exit(rm("control", envir=env))
  UseMethod("nmObjHandleControlObject")
}

#' @rdname nmObjHandleControlObject
#' @export
nmObjHandleControlObject.foceiControl <- function(control, env) {
  assign("foceiControl", control, envir=env)
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
  if (exists("foceiControl", .env)) {
    .control <- get("foceiControl", .env)
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
nmObjGetControl.foce <- nmObjGetControl.focei

#' @rdname nmObjGetControl
#' @export
nmObjGetControl.foi <- nmObjGetControl.focei

#' @rdname nmObjGetControl
#' @export
nmObjGetControl.fo <- nmObjGetControl.focei

#' @rdname nmObjGetControl
#' @export
nmObjGetControl.posthoc <- nmObjGetControl.focei

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
