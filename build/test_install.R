                                        # Installation test function
test_install <- function(){
                                        # Test 1: Correct R version
    if(sessionInfo()$R.version$major=="3" & as.numeric(sessionInfo()$R.version$minor)>=4.1){
        cat(paste0("Correct R version: Yes, ",sessionInfo()$R.version$version.string,"\n"))
    }else{
        cat("Correct R version: No, go to https://cloud.r-project.org/ and install latest version\n")
        return("Installation not complete!")
    }
    tmp <- try(library(rxode2), silent=TRUE);
    if (inherits(tmp, "try-error"))  {
        cat("rxode2 installed: No; check https://github.com/nlmixr2/nlmixr2 for instructions to do so\n")
        return("Installation not complete!")
    } else {
        cat(paste0("rxode2 installed: Yes\n"))
    }

                                        # Test 4: devtools installed
    if("devtools" %in% rownames(installed.packages())){
        cat("devtools package installed: Yes\n")
    }else{
        cat("devtools package installed: No, use install.packages('devtools') to do so\n")
        return("Installation not complete!")
    }

                                        # Test 6: rxode2 installed
    if("rxode2" %in% rownames(installed.packages())){
        cat("rxode2 package installed: Yes\n")
    }else{
        cat("rxode2 package installed: No, use install_github('nlmixr2/rxode2') to do so\n")
        return("Installation not complete!")
    }

    ## Test 5: Rtools installed
    if(rxode2:::.rxWinRtoolsPath()){
        cat("Rtools installed or not needed: Yes\n")
    }else{
        cat("Rtools installed: No, check https://github.com/nlmixr2/nlmixr2 for instructions to do so\n")
        return("Installation not complete!")
    }

    if ("units" %in% rownames(installed.packages())){
        cat("units package installed: Yes\n")
    } else {
        cat("units package not installed.\n");
        cat("Units requires udunits on some systems.\n\tUse https://github.com/r-quantities/units#installation for installation instructions.\n");
        return("Installation not complete!")
    }

                                        # Test 8: nlmixr installed
    if("nlmixr" %in% rownames(installed.packages())){
        cat("nlmixr package installed: Yes\n")
    }else{
        cat("nlmixr package installed: No, use install_github('nlmixr2/nlmixr2') to do so\n")
        return("Installation not complete!")
    }

                                        # Test 9: xpose.nlmixr installed
    if("xpose.nlmixr" %in% rownames(installed.packages())){
        cat("xpose.nlmixr package installed: Yes\n")
    }else{
        cat("xpose.nlmixr package installed: No, use install_github('nlmixr2/xpose.nlmixr') to do so\n")
        return("Installation not complete!")
    }

                                        # Test 10: shinyMixR installed
    if("shinyMixR" %in% rownames(installed.packages())){
        cat("shinyMixR package installed: Yes\n")
    }else{
        cat("shinyMixR package installed: No, use install_github('richardhooijmaijers/shinyMixR') to do so\n")
        return("Installation not complete!")
    }

                                        # Test 11: nlme test for theophylline
    library(nlmixr)
    testmod <- function() {
      ini({
        tka <- 0.45 # Log Ka
        tcl <- 1 # Log Cl
        ## This works with interactive models
        ## You may also label the preceding line with label("label text")
        tv <- 3.45; label("log V")
        ## the label("Label name") works with all models
        eta.ka ~ 0.6
        eta.cl ~ 0.3
        eta.v ~ 0.1
        add.sd <- 0.7
      })
      model({
        ka <- exp(tka + eta.ka)
        cl <- exp(tcl + eta.cl)
        v <- exp(tv + eta.v)
        linCmt() ~ add(add.sd)
      })
    }
    if(Sys.info()['sysname']%in%c('Darwin','Linux')){
        message("Running nlme")
        {sink("/dev/null"); fit1 <- suppressWarnings(suppressMessages(nlmixr(testmod,theo_sd,est="nlme"))); sink(); }
        message("Running saem")
        {sink("/dev/null"); fit2 <- suppressWarnings(suppressMessages(nlmixr(testmod,theo_sd,est="saem"))); sink(); }
        message("Running focei")
        {sink("/dev/null"); fit3 <- suppressWarnings(suppressMessages(nlmixr(testmod,theo_sd,est="focei"))); sink(); }
    }else{
        message("Running nlme")
        {sink("NUL"); fit1 <- suppressWarnings(suppressMessages(nlmixr(testmod,theo_sd,est="nlme"))); sink(); }
        message("Running saem")
        {sink("NUL"); fit2 <- suppressWarnings(suppressMessages(nlmixr(testmod,theo_sd,est="saem"))); sink(); }
        message("Running focei")
        {sink("NUL"); fit3 <- suppressWarnings(suppressMessages(nlmixr(testmod,theo_sd,est="focei"))); sink(); }
    }
    if(is(fit1, "nlmixrFitData") && !is.null(fit1$nlme)){
        cat("nlmixr run under nlme: Yes\n")
    }else{
        cat("nlmixr run under nlme: No, contact nlmixr team\n")
        return("Installation not complete!")
    }
    if(is(fit2, "nlmixrFitData") && !is.null(fit2$saem)){
        cat("nlmixr run under saem: Yes\n")
    }else{
        cat("nlmixr run under saem: No, contact nlmixr team\n")
        return("Installation not complete!")
    }
    if(is(fit3, "nlmixrFitData") && fit3$method == "FOCE"){
        cat("nlmixr run under focei: Yes\n")
    }else{
        cat("nlmixr run under focei: No, contact nlmixr team\n")
        return("Installation not complete!")
    }
    cat("---- Installation test finished! ----\n")
}

try({test_install()})

cat("Begin Session Info:\n")
print(sessionInfo())
cat("\n")
