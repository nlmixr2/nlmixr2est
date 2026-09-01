## Regenerate tests/testthat/baselines/qrpem-baseline-ref.rds.
##
## The baseline pins that impmap's UN-ADAPTED path (auto=FALSE, the one-sided
## gammaRule) is unchanged, so its value is in being old -- do NOT regenerate it
## to make a failing diff go away.  Regenerate ONLY when a change is understood
## to move this fit on purpose, and record which commit did it in the test.
##
##   Rscript inst/tools/genQrpemBaseline.R
##
## Run from the package root against the build you want recorded.
suppressMessages(devtools::load_all(".", helpers = FALSE, quiet = TRUE))

.oneCmt <- function() {
  ini({
    tka <- 0.45; tcl <- 1; tv <- 3.45
    eta.ka ~ 0.6; eta.cl ~ 0.3
    add.sd <- 0.7
  })
  model({
    ka <- exp(tka + eta.ka)
    cl <- exp(tcl + eta.cl)
    v <- exp(tv)
    linCmt() ~ add(add.sd)
  })
}

## MUST match the control in test-qrpem-slow.R exactly.
.f <- suppressWarnings(
  nlmixr2(.oneCmt, nlmixr2data::theo_sd, "impmap",
          impmapControl(print = 0L, nIter = 5L, isample = 100L, auto = FALSE,
                        sigdig = 4, gammaRule = "floor")))

saveRDS(list(fixef = fixef(.f),
             omega = .f$omega,
             obj = .f$env$impObj,
             samples1 = .f$env$impSamples[[1]]),
        file.path("tests", "testthat", "baselines", "qrpem-baseline-ref.rds"))
cat("wrote tests/testthat/baselines/qrpem-baseline-ref.rds\n")
print(fixef(.f))
