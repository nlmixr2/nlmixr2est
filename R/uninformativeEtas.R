#' Get the parameter values for
#'
#' @param name name of the eta parameter
#' @param ui rxode2 user interface object
#' @param pm plus or minus values
#' @param plus What type of deviate is this:
#'
#' - `plus=TRUE`: high estimate
#' - `plus=FALSE`: low estimate
#' - `plus=NA`: middle estimate
#' @param saem is this a saem-style mu-referenced model?
#' @return the named parameter value for plus/minus/mid calculation
#' @noRd
#' @author Matthew L. Fidler
.getMuValForUE <- function(name, ui, pm, plus=TRUE, saem=TRUE, retName=FALSE) {
  if (!saem || name %in% ui$nonMuEtas) {
    if (retName) return(name)
    if (is.na(plus)) {
      return(setName(0, name))
    } else if (plus) {
      return(pm[name])
    } else {
      return(-pm[name])
    }
  }
  .w <- which(ui$muRefDataFrame$eta == name)
  .n2 <- ui$muRefDataFrame$theta[.w]
  .v0 <- ui$theta[.n2]
  if (retName) return(.n2)
  if (is.na(plus)) {
    setNames(.v0, .n2)
  } else if (plus) {
    setNames(.v0 + pm[name], .n2)
  } else {
    setNames(.v0 - pm[name], .n2)
  }
}
#' Get the low, middle, and hi subject values based on etas
#'
#' @param ui user interface function
#' @param data data.frame that will be used for fitting
#' @param alpha The alpha value to scale the interval from -1 to 1 when choosing the quadrature points
#' @param saem boolean indicating if the model is a mu-referenced.
#' @param q plus or minus quadrature
#' @return a list with:
#' - `low` lower parameters
#' - `mid` middle individual parameters
#' - `hi` high individual parameters
#' @noRd
#' @author Matthew L. Fidler
.uninformativeEtas <- function(ui, data, alpha=0.05,
                               saem=TRUE, q=sqrt(3/5)) {
  ui <- rxode2::assertRxUi(ui)
  .trans <- rxode2::etTrans(data, ui)
  .lst <- attr(class(.trans), ".rxode2.lst")
  .n <- .lst$nid
  .iniDf <- ui$iniDf
  .eta <- .iniDf[which(.iniDf$neta1 == .iniDf$neta2), ]

  .pm <- setNames(qnorm(1 - alpha / 2) * sqrt(.eta$est) * q, .eta$name)

  .nn <- vapply(names(.pm), .getMuValForUE, ui=ui, pm=.pm, plus=TRUE, saem=saem, retName=TRUE,
                character(1), USE.NAMES=FALSE)

  .p <- do.call("rbind", lapply(seq_len(.n), function(i) {
    as.data.frame(t(vapply(names(.pm), .getMuValForUE, ui=ui, pm=.pm, plus=TRUE, saem=saem,
                           double(1), USE.NAMES=FALSE)))
  }))
  names(.p) <- .nn

  .m <- do.call("rbind", lapply(seq_len(.n), function(i) {
    as.data.frame(t(vapply(names(.pm), .getMuValForUE, ui=ui, pm=.pm, plus=FALSE, saem=saem,
                           double(1), USE.NAMES = FALSE)))
  }))
  names(.m) <- .nn

  .z <- do.call("rbind", lapply(seq_len(.n), function(i) {
    as.data.frame(t(vapply(names(.pm), .getMuValForUE, ui=ui, pm=.pm, plus=NA, saem=saem,
                           double(1), USE.NAMES = FALSE)))
  }))
  names(.z) <- .nn

  list(trans=setNames(names(.pm), .nn), low=.m, mid=.z, hi=.p)
}
