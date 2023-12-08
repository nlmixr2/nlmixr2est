## nmTest({

##   rxode2::rxUnloadAll()
##   gc()

##   dat <- Wang2007
##   dat$DV <- dat$Y

##   gg <- function(x, y) {
##     x * y
##   }

##   f <- function() {
##     ini({
##       tke <- 0.5
##       eta.ke ~ 0.04
##       prop.sd <- sqrt(0.1)
##     })
##     model({
##       ke <- gg(tke, exp(eta.ke))
##       ipre <- gg(10, exp(-ke * t))
##       lipre <- log(ipre)
##       ipre ~ prop(prop.sd)
##     })
##   }

##   .env <- new.env(parent=emptyenv())

##   expect_error({.env$f <- nlmixr2(f, dat, "focei")}, NA)

##   udfCheck <- function(udf) {
##     .w <- which(is.na(udf))
##     udf[-.w]
##   }

##   expect_equal(udfCheck(rxModelVars(.env$f$foceiModel$inner)$udf), c("gg"=2L))
##   expect_equal(udfCheck(rxModelVars(.env$f$foceiModel$predOnly)$udf), c("gg"=2L))
##   expect_equal(udfCheck(rxModelVars(.env$f$foceiModel$predNoLhs)$udf), c("gg"=2L))

##   rxode2::rxUnloadAll()
##   gc()

##   expect_error({.env$s <- nlmixr2(f, dat, "saem")}, NA)

##   rxode2::rxUnloadAll()
##   gc()

##   expect_error({.env$n <- nlmixr2(f, dat, "nlme")}, NA)

##   g <- function() {
##     ini({
##       tke <- 0.5
##       add.sd <- sqrt(0.1)
##     })
##     model({
##       ke <- tke
##       ipre <- gg(10, exp(-ke * t))
##       lipre <- log(ipre)
##       ipre ~ add(add.sd)
##     })
##   }

##   rxode2::rxUnloadAll()
##   gc()

##   expect_error({.env$nlm <- nlmixr2(g, dat, "nlm")}, NA)

##   rxode2::rxUnloadAll()
##   gc()


##   expect_error({.env$optim <- nlmixr2(g, dat, "optim")}, NA)

##   rxode2::rxUnloadAll()
##   gc()


##   ## expect_error({.env$nls <- nlmixr2(g, dat, "nls")}, NA)

##   expect_error({.env$nlminb <- nlmixr2(g, dat, "nlminb")}, NA)

##   rxode2::rxUnloadAll()
##   gc()


##   expect_error({.env$bobyqa <- nlmixr2(g, dat, "bobyqa")}, NA)

##   rxode2::rxUnloadAll()
##   gc()

##   expect_error({.env$lbfgsb3c <- nlmixr2(g, dat, "lbfgsb3c")}, NA)

##   rxode2::rxUnloadAll()
##   gc()

##   expect_error({.env$n1qn1 <- nlmixr2(g, dat, "n1qn1")}, NA)

##   rxode2::rxFun(gg)

##   rm(gg)

##   rxClean()

##   rxode2::rxUnloadAll()
##   gc()

##   .env$f2 <- nlmixr2(f, dat, "focei")

##   rxode2::rxUnloadAll()
##   gc()

##   .env$s2 <- nlmixr2(f, dat, "saem")

##   rxode2::rxUnloadAll()
##   gc()

##   .env$n2 <- nlmixr2(f, dat, "nlme")

##   rxode2::rxUnloadAll()
##   gc()


##   rxode2::rxRmFun("gg")

## })
