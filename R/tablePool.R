#' Set up the post-fit table models in the ODE-swap pool
#'
#' Builds one solve pool over the fit's translated events, sized for the larger
#' of `inner` and `predOnly`, and registers both as pool slots (#1356).
#'
#' @param fit focei style fit
#' @param inner the table's ipred/pred model (`fit$innerModel`), or `NULL`
#' @param predOnly the display-column model (`fit$predOnlyModel`), or `NULL`
#' @param pars parameters the pool is built with (one row per subject)
#' @param method rxode2 ODE method name
#' @param atol,rtol ODE tolerances
#' @return named logical: whether each model can be solved in the pool
#' @author Matthew L. Fidler
#' @noRd
.tablePoolSetup <- function(fit, inner, predOnly, pars, method, atol, rtol) {
  .pool <- tablePoolDeclare_(inner, predOnly)
  if (is.null(.pool)) {
    return(c(inner = FALSE, predOnly = FALSE))
  }
  .covsi <- .residCovsInterpolation(fit)
  .rxSolveNonmem(
    .pool,
    pars,
    fit$dataSav,
    atol = atol,
    rtol = rtol,
    maxsteps = fit$maxstepsOde,
    hmin = fit$hmin,
    hmax = fit$hmax,
    hini = fit$hini,
    maxordn = fit$maxordn,
    maxords = fit$maxords,
    method = rxode2::odeMethodToInt(method),
    tolFactor = fit$env$tolFactor,
    covsInterpolation = .covsi$covsInterpolation,
    naInterpolation = .covsi$naInterpolation,
    nonmem = .covsi$nonmem,
    warnIdSort = FALSE,
    .setupOnly = 1L
  )
  tablePoolRegister_(inner, predOnly)
}

#' Solve one table model in the pool
#'
#' @param model the model, which must be the one registered for `slot`
#' @param slot `"inner"` or `"predOnly"`
#' @param pars named parameter vector (all subjects) or matrix (one row per
#'   subject)
#' @param lhs,state lhs/state names to output as columns
#' @param sens lhs names to return as sparse triplets instead
#' @param cov model parameters (covariates) to output, as interpolated
#' @param tbs output `rxLambda`/`rxYj`/`rxLow`/`rxHi`
#' @return list: `id` (subject index), `rec` (0-based record within the
#'   subject), `evid`, `time`, `vals` (data frame of the requested columns),
#'   `sens` (triplets) and `rowsOk`
#' @author Matthew L. Fidler
#' @noRd
.tablePoolSolve <- function(
  model,
  slot,
  pars,
  lhs = character(0),
  state = character(0),
  sens = character(0),
  cov = character(0),
  tbs = FALSE
) {
  .mv <- rxode2::rxModelVars(model)
  if (is.null(dim(pars))) {
    pars <- matrix(pars, nrow = 1L, dimnames = list(NULL, names(pars)))
  } else {
    pars <- as.matrix(pars)
  }
  .idx <- function(x, from) {
    .i <- match(x, from)
    if (anyNA(.i)) {
      stop("not in the table model: ", paste(x[is.na(.i)], collapse = ", "), call. = FALSE)
    }
    as.integer(.i - 1L)
  }
  .ret <- tablePoolSolve_(
    switch(slot, inner = 0L, predOnly = 1L),
    pars,
    .idx(colnames(pars), .mv$params),
    .idx(lhs, .mv$lhs),
    .idx(sens, .mv$lhs),
    .idx(state, .mv$state),
    .idx(cov, .mv$params),
    tbs
  )
  .vals <- as.data.frame(.ret$vals)
  names(.vals) <- c(lhs, state, cov, if (tbs) c("rxLambda", "rxYj", "rxLow", "rxHi"))
  .ret$vals <- .vals
  .ret
}
