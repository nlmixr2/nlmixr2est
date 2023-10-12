#' nlmixr2 defaults controls for nlm
#'
#' @inheritParams stats::nlm
#' @inheritParams nlmixrNlmeControl
#' @return
#' @export
#' @author Matthew L. Fidler
#' @examples
nlmControl <- function(hessian = TRUE, typsize = NULL,
                       fscale = 1, print.level = 2, ndigit = NULL, gradtol = 1e-6,
                       stepmax = NULL,
                       steptol = 1e-6, iterlim = 10000, check.analyticals = TRUE,
                       rxControl=NULL,
                       optExpression=TRUE, sumProd=FALSE,
                       returnNlm=FALSE,
                       addProp = c("combined2", "combined1"),
                       calcTables=TRUE, compress=TRUE,
                       adjObf=TRUE, ci=0.95, sigdig=4, sigdigTable=NULL, ...) {
  checkmate::assertLogical(hessian, len=1, any.missing=FALSE)
  checkmate::assertNumeric(stepmax, lower=0, len=1, null.ok=TRUE, any.missing=FALSE)
  checkmate::assertIntegerish(print.level, lower=1, upper=2, any.missing=FALSE)
  checkmate::assertNumeric(ndigit, lower=0, len=1, any.missing=FALSE, null.ok=TRUE)
  checkmate::assertNumeric(gradtol, lower=0, len=1, any.missing=FALSE)
  checkmate::assertNumeric(steptol, lower=0, len=1, any.missing=FALSE)
  checkmate::assertIntegerish(iterlim, lower=1, len=1, any.missing=FALSE)
  checkmate::assertLogical(check.analyticals, len=1, any.missing=FALSE)
  checkmate::assertLogical(returnNlm, len=1, any.missing=FALSE)
  checkmate::assertLogical(calcTables, len=1, any.missing=FALSE)
  checkmate::assertLogical(compress, len=1, any.missing=TRUE)
  checkmate::assertLogical(adjObf, len=1, any.missing=TRUE)
  .xtra <- list(...)
  .bad <- names(.xtra)
  .bad <- .bad[!(.bad %in% c("genRxControl"))]
  if (length(.bad) > 0) {
    stop("unused argument: ", paste
    (paste0("'", .bad, "'", sep=""), collapse=", "),
    call.=FALSE)
  }

  .genRxControl <- FALSE
  if (!is.null(.xtra$genRxControl)) {
    .genRxControl <- .xtra$genRxControl
  }
  if (is.null(ndigit)) {
    ndigit <- sigdig
  }
  if (is.null(rxControl)) {
    if (!is.null(sigdig)) {
      rxControl <- rxode2::rxControl(sigdig=sigdig)
    } else {
      rxControl <- rxode2::rxControl(atol=1e-4, rtol=1e-4)
    }
    .genRxControl <- TRUE
  } else if (inherits(rxControl, "rxControl")) {
  } else if (is.list(rxControl)) {
    rxControl <- do.call(rxode2::rxControl, rxControl)
  } else {
    stop("solving options 'rxControl' needs to be generated from 'rxode2::rxControl'", call=FALSE)
  }
  if (!is.null(sigdig)) {
    checkmate::assertNumeric(sigdig, lower=1, finite=TRUE, any.missing=TRUE, len=1)
    if (is.null(sigdigTable)) {
      sigdigTable <- round(sigdig)
    }
  }
  if (is.null(sigdigTable)) {
    sigdigTable <- 3
  }
  checkmate::assertIntegerish(sigdigTable, lower=1, len=1, any.missing=FALSE)

  .ret <- list(hessian = hessian, typsize = typsize,
               fscale = fscale, print.level = print.level, ndigit=ndigit, gradtol = gradtol,
               stepmax = stepmax,
               steptol = steptol, iterlim = iterlim,
               check.analyticals = check.analyticals,
               optExpression=optExpression,
               sumProd=sumProd,
               rxControl=rxControl,
               returnNlm=returnNlm, addProp=addProp, calcTables=calcTables,
               compress=compress,
               ci=ci, sigdig=sigdig, sigdigTable=sigdigTable,
               genRxControl=.genRxControl)
  class(.ret) <- "nlmControl"
  .ret
}
