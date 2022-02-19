#' Assert that this is a nlmixr2 fit object
#'
#' Will error without nlmixr2 fit object
#'
#' @param fit Fit object
#' @return Nothing
#' @author Matthew L. Fidler
#' @examples
#' \dontrun{
#'
#' f <- 4
#' assertNlmixrFit(f) # throw error
#'
#' }
#' @export
assertNlmixrFit <- function(fit) {
  .object <- as.character(substitute(fit))
  if (!(inherits(fit, "nlmixr2FitCore") ||
        inherits(fit, "nlmixr2FitCoreSilent"))) {
    stop("'", .object, "' needs to be a nlmixr2 fit object",
         call.=FALSE)
  }
}

#' Assert that this is a nlmixr2 fit data object
#'
#' Will error without nlmixr2 fit data object
#'
#' @param fit Fit object
#' @return Nothing
#' @author Matthew L. Fidler
#' @examples
#' \dontrun{
#'
#' f <- 4
#' assertNlmixrFitData(f) # throw errors
#'
#' }
#' @export
assertNlmixrFitData <- function(fit) {
  .object <- as.character(substitute(fit))
  if (!inherits(fit, "nlmixr2FitData")) {
    stop("'", .object, "' needs to be a nlmixr2 fit object with data attached",
         call.=FALSE)
  }
}
#' Assert a nlmixr2 object data frame row is compatible with what needs to be added
#'
#' @param df Data frame to assert
#' @param allowNa Allow NA data frame
#' @return If successful, a list(data frame, condition number)
#' @author Matthew L. Fidler
#' @noRd
assertNlmixrObjDataFrameRow <- function(df, allowNa=FALSE) {
  .name <- names(df)
  .needed <- c("OBJF", "AIC", "BIC", "Log-likelihood")
  .diff <- setdiff(.needed, .name)
  if (length(.diff) > 0) {
    stop("need additional information for objective function data frame row: '", paste(.diff, collapse="', '"), "'",
         call.=FALSE)
  }
  .w <- which(.name == "Condition Number")
  if (length(.w) == 1) {
    .cn <- df[["Condition Number"]]
    if (inherits(.cn, "numeric")) {
      if (!is.na(.cn)) {
        .cn <- setNames(.cn, NULL)
      } else {
        .cn <- NA
      }
    } else {
      .cn <- NA
    }
  } else {
    .cn <- NA
  }
  .df1 <- df[, .needed]
  if (is.na(.df1[, 1])) {
    .df1 <- NA
    if (!allowNa) stop("missing objective function in objective function data frame", call.=FALSE)
  } else {
    lapply(.needed, function(x) {
      checkmate::assertNumeric(df[[x]], len=1, any.missing=FALSE, .var.name=x)
    })
  }
  list(.df1, .cn)
}
