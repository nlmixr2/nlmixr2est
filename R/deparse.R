.deparseShared <- function(x, value) {
  if (x == "rxControl") {
    .rx <- rxUiDeparse(value, "a")
    .rx <- .rx[[3]]
    return(paste0("rxControl = ", deparse1(.rx)))
  } else if (x == "scaleType")  {
    .scaleTypeIdx <- c("norm" = 1L, "nlmixr2" = 2L, "mult" = 3L, "multAdd" = 4L)
    paste0("scaleType =", deparse1(names(.scaleTypeIdx[which(object[[x]] == .scaleTypeIdx)])))
  } else if (x == "normType") {
    .normTypeIdx <- c("rescale2" = 1L, "rescale" = 2L, "mean" = 3L, "std" = 4L, "len" = 5L, "constant" = 6L)
    paste0("normType =", deparse1(names(.normTypeIdx[which(object[[x]] == .normTypeIdx)])))
  }
  return(NA_character_)
}

.deparseDifferent <- function(standard, new, internal=character(0)) {
  which(vapply(names(standard),
               function(x) {
                 if (is.function(standard[[x]])) {
                   warning(paste0(x, "as a function not supported in ",
                                  class(standard), "()"), call.=FALSE)
                   FALSE
                 } else if (x %in% internal){
                   FALSE
                 } else {
                   !identical(standard[[x]], new[[x]])
                 }
               }, logical(1), USE.NAMES=FALSE))
}

.deparseFinal <- function(default, object, w, var, fun=NULL) {
  .cls <- class(object)
  if (length(w) == 0) {
    return(str2lang(paste0(var, " <- ", .cls, "()")))
  }
  .retD <- vapply(names(default)[w], function(x) {
    .val <- .deparseShared(x, object[[x]])
    if (!is.na(.val)) {
      return(.val)
    }
    if (is.function(fun)) {
      .val <- fun(default, x, object[[x]])
      if (!is.na(.val)) {
        return(.val)
      }
    }
    paste0(x, "=", deparse1(object[[x]]))
  }, character(1), USE.NAMES=FALSE)
  str2lang(paste(var, " <- ", .cls, "(", paste(.retD, collapse=","),")"))
}
