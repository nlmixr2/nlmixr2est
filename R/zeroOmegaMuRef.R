# Mu-referenced random effects whose variance is fixed at zero.
#
# This is a computational device, not a model feature.  A control stream writes
# `MU_5 = THETA(5)` with `$OMEGA (0.0 FIXED)` to give THETA(5) an EM update
# path -- the random effect carries no variability of its own, it exists so the
# method has something to take a conditional mean of.  Bauer's non-Gaussian-eta
# streams do this for all five distribution parameters.
#
# Pinning such a variance at zero the way an ordinary fix() is pinned does not
# work, and fails silently:
#
#   * saem's mu-theta M-step is a normal equation weighted by omega^-1
#     (`Plambda1`, src/saem.cpp), so a zero variance gives that column infinite
#     weight, the sampler cannot move it off its prior mean -- which IS the
#     current theta -- and the M-step's fixed point becomes the theta's
#     STARTING value.  Measured on Bauer's gamma model: CL stayed at 6.686
#     against a truth of 5.03, reported as an estimate.
#   * imp's E-step gets its proposal from the MAP, and a zero prior variance
#     puts the MAP at eta = 0 for every subject, so the mean conditional eta
#     that impMuInterceptStep() folds into the theta is identically zero and
#     the theta never moves either.
#
# Both methods already have the right machinery -- the EM update for such a
# theta is "shift it by the mean of its random effect's conditional values"
# (impMuInterceptStep() in src/inner.cpp does exactly that, and saem's
# regression does it implicitly).  What that machinery needs is for the random
# effect to be ABLE to move, which means letting the variance estimate while
# the method runs and treating it as the working quantity it is.  The variance
# is then reported as the zero it was declared to be, because by then the
# theta has absorbed everything it was carrying.
#
# @param ui rxode2 model
# @return character vector of eta names to treat this way (possibly empty)
# @noRd
# @author Matthew L. Fidler
.zeroOmegaMuRefEtas <- function(ui) {
  .iniDf <- ui$iniDf
  if (is.null(.iniDf) || !any(names(.iniDf) == "neta1")) return(character(0))
  .eta <- .iniDf[!is.na(.iniDf$neta1), ]
  if (nrow(.eta) == 0L) return(character(0))
  .diag <- .eta[.eta$neta1 == .eta$neta2, ]
  ## Fixed at (numerically) zero.  Below this the normal equations are
  ## conditioned past what a double carries, so the pin is degenerate whatever
  ## the user meant by the value.
  .cand <- .diag[.diag$fix & abs(.diag$est) <= 1e-8, ]
  if (nrow(.cand) == 0L) return(character(0))
  ## A zero variance cannot carry a covariance, so anything with an
  ## off-diagonal is not this construct and is left alone.
  .off <- .eta[.eta$neta1 != .eta$neta2, ]
  if (nrow(.off) > 0L) {
    .bad <- unique(c(.off$neta1, .off$neta2))
    .cand <- .cand[!(.cand$neta1 %in% .bad), ]
  }
  if (nrow(.cand) == 0L) return(character(0))
  ## Only mu-referenced ones: the whole point is that a theta absorbs the
  ## random effect's location, and without a mu reference there is no theta to
  ## absorb it into.
  .mr <- ui$muRefDataFrame
  if (is.null(.mr) || nrow(.mr) == 0L) return(character(0))
  intersect(.cand$name, .mr$eta)
}

## Methods whose theta update for a mu-referenced parameter is "shift by the
## mean of its random effect".  Only these have anything to gain here.  The
## FOCEi family estimates such a theta directly through its outer optimizer, so
## giving the random effect a working variance there would just add
## between-subject variability the model does not have.
.zeroOmegaMuRefMethods <- c("saem", "imp", "impmap", "npag", "npb")

#' Read back the etas the preprocess hook rewrote
#' @param ui rxode2 model (the rewritten one carried by the fit)
#' @return character vector, possibly empty
#' @noRd
.zeroOmegaMuRefStash <- function(ui) {
  if (is.null(ui) || !is.environment(ui$meta)) return(character(0))
  if (!exists(".zeroOmegaMuRefEtas", envir=ui$meta)) return(character(0))
  get(".zeroOmegaMuRefEtas", envir=ui$meta)
}

#' The thetas those random effects are mu-referenced to
#' @param ui rxode2 model carried by the fit
#' @return character vector, possibly empty
#' @noRd
.zeroOmegaMuRefThetaStash <- function(ui) {
  if (is.null(ui) || !is.environment(ui$meta)) return(character(0))
  if (!exists(".zeroOmegaMuRefThetas", envir=ui$meta)) return(character(0))
  get(".zeroOmegaMuRefThetas", envir=ui$meta)
}

#' Let a zero-variance mu-referenced random effect move while the method runs
#'
#' See .zeroOmegaMuRefEtas() for why the declared zero cannot simply be held.
#' The variance is seeded with a working value and un-fixed, so it is an
#' ordinary estimated random effect for the duration of the fit; every method's
#' existing machinery then does the right thing with it, and the declared zero
#' is put back by .postFinalZeroOmegaMuRef() once the theta has absorbed what it
#' was carrying.
#'
#' @param ui rxode2 model
#' @param est estimation method name
#' @param data input data
#' @param control method control
#' @return list with the rewritten `ui`, or NULL when there is nothing to do
#' @noRd
#' @author Matthew L. Fidler
.preProcessZeroOmegaMuRef <- function(ui, est, data, control) {
  if (!is.character(est) || length(est) != 1L) return(NULL)
  if (!(est %in% .zeroOmegaMuRefMethods)) return(NULL)
  .zero <- .zeroOmegaMuRefEtas(ui)
  if (length(.zero) == 0L) return(NULL)
  .ui <- rxode2::rxUiDecompress(ui)
  .iniDf <- .ui$iniDf
  .w <- which(.iniDf$name %in% .zero &
                !is.na(.iniDf$neta1) & .iniDf$neta1 == .iniDf$neta2)
  if (length(.w) == 0L) return(NULL)
  ## What the two families need here is NOT the same, so this is split by
  ## method rather than done one way for both.
  ##
  ## The FOCEi family (focei/imp/impmap/npag) must see the declared ZERO: its
  ## omega setup detects a flat random effect by looking for a zero diagonal
  ## and substitutes -- then removes -- its own placeholder.  Rewriting the
  ## zero here hides it from exactly that code (measured: `fix(0)` came through
  ## with no flat indices while `fix(1e-9)`, the same model, was excluded).
  ##
  ## saem has no such local step: its omega is consumed through covstruct,
  ## Gamma2_phi1 and the MCMC before any of that, and a zero there means the
  ## random effect cannot move at all -- so the conditional mean the M-step
  ## folds into the theta is identically zero and the theta never budges.  It
  ## gets the exploration value written in, held fix()ed so the sufficient
  ## statistics do not update it, and reported back as the declared zero.
  ## saem's omega IS rewritten here, and the value comes from the `control`
  ## ARGUMENT rather than from the model.
  ##
  ## This hook runs before the method's control is attached to the ui, so
  ## reading it with rxGetControl(ui, ...) silently returns the default -- and
  ## that silently disabled `saemControl(zeroOmegaTune=)` entirely: the hook
  ## wrote 0.1 whatever the user asked for, `.configsaem()`'s own substitution
  ## then found no zeros left to replace, and the value it was correctly handed
  ## did nothing.  Measured: tune = 0.1 and tune = 4 gave bit-identical fits.
  ##
  ## The rewrite cannot simply be dropped in favour of the one in
  ## `.configsaem()`: without it saem returns Inf with the theta stuck at its
  ## starting value.  Its omega is consumed through covstruct and the MCMC
  ## before .configsaem's substitution can matter.
  if (identical(est, "saem")) {
    .tune <- NULL
    if (is.list(control) && !is.null(control$zeroOmegaTune)) {
      .tune <- control$zeroOmegaTune
    }
    if (is.null(.tune)) {
      .tune <- rxode2::rxGetControl(.ui, "zeroOmegaTune", 0.1)
    }
    if (!is.numeric(.tune) || length(.tune) != 1L || !is.finite(.tune) ||
          .tune <= 0) {
      .tune <- 0.1
    }
    .iniDf$est[.w] <- .tune
    .iniDf$fix[.w] <- TRUE
    assign("iniDf", .iniDf, envir=.ui)
  }
  ## Recorded rather than recomputed later: after this rewrite the model no
  ## longer looks like the thing that triggered it -- the variance is an
  ## ordinary estimated one and there is nothing left to detect.  It goes in
  ## the model's meta environment because that is what survives to the fit;
  ## rxAssignControlValue() does not (measured: readable on the ui the hook
  ## returns, gone by the time the fit object exists), which silently left the
  ## working variance reported as an estimate.
  ## The THETA each of these random effects is mu-referenced to is recorded
  ## here, not re-derived at fit time.  A zero-variance random effect is
  ## missing from `muRefDataFrame` by the time a method looks (measured: with
  ## `fix(0)` the later lookup matched nothing and the flat columns never
  ## reached saem at all, while `fix(1e-9)` -- the same model under this
  ## rewrite -- came through), so the pairing has to be taken while the
  ## original model is still in hand.
  .mr <- ui$muRefDataFrame
  .th <- character(0)
  if (!is.null(.mr) && nrow(.mr) > 0L) {
    .th <- .mr$theta[.mr$eta %in% .zero]
  }
  if (is.environment(.ui$meta)) {
    assign(".zeroOmegaMuRefEtas", .zero, envir=.ui$meta)
    assign(".zeroOmegaMuRefThetas", .th, envir=.ui$meta)
  }
  list(ui=.ui)
}

preProcessHooksAdd(".preProcessZeroOmegaMuRef", .preProcessZeroOmegaMuRef)

#' Report a zero-variance mu-referenced random effect as the zero it is
#'
#' The variance was left estimable while the method ran so the random effect
#' could move and the EM had a conditional mean to shift the theta by.  By the
#' time the fit is assembled the theta has absorbed what it was carrying, so
#' the working variance is put back to the declared zero -- reporting it as an
#' estimate would claim between-subject variability the model does not have.
#'
#' `saem` additionally does this in `.getSaemOmega()`, which builds its own
#' `$omega` while keeping saem's phi ordering straight.
#'
#' @param ret fit object
#' @return `ret`, with the declared zeros restored
#' @noRd
#' @author Matthew L. Fidler
.postFinalZeroOmegaMuRef <- function(ret) {
  .ui <- try(ret$ui, silent=TRUE)
  if (inherits(.ui, "try-error") || is.null(.ui)) return(ret)
  .zero <- try(.zeroOmegaMuRefStash(.ui), silent=TRUE)
  if (inherits(.zero, "try-error") || length(.zero) == 0L) return(ret)
  .om <- try(ret$omega, silent=TRUE)
  if (inherits(.om, "try-error") || is.null(.om) || is.null(dimnames(.om))) {
    return(ret)
  }
  .w <- which(rownames(.om) %in% .zero)
  if (length(.w) == 0L) return(ret)
  for (.k in .w) {
    .om[.k, ] <- 0
    .om[, .k] <- 0
  }
  assign("omega", .om, envir=ret$env)
  ret
}

postFinalObjectHooksAdd(".postFinalZeroOmegaMuRef", .postFinalZeroOmegaMuRef)
