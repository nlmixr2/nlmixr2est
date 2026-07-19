nmTest({
  # Mu-referenced thetas (pop + covariate coefficients) are regression-updated
  # (no outer-optimizer slot) but appear as standard scale.h columns at their
  # natural theta positions and in parHist/parHistData; gradient rows record
  # NaN (blank console cells). Console layout details are covered in
  # test-mfocei.R (weekly batch); this is the always-run core check.

  .ocmt <- function() {
    ini({
      tka <- 0.45
      tcl <- 1
      tv <- 3.45
      eta.ka ~ 0.6
      eta.cl ~ 0.3
      eta.v ~ 0.1
      add.sd <- 0.7
    })
    model({
      ka <- exp(tka + eta.ka)
      cl <- exp(tcl + eta.cl)
      v <- exp(tv + eta.v)
      d/dt(depot) <- -ka * depot
      d/dt(center) <- ka * depot - cl / v * center
      cp <- center / v
      cp ~ add(add.sd)
    })
  }

  test_that("ifocei parHist has mu-theta columns in natural order with NaN gradients", {
    skip_on_cran()

    d <- nlmixr2data::theo_sd
    out <- withr::with_options(list(width = 200), capture.output({
      fitM <- nlmixr2est::nlmixr(.ocmt, d, "ifocei",
                                 ifoceiControl(print = 1, maxOuterIterations = 2L,
                                               outerOpt = "nlminb", covMethod = "", calcTables = FALSE))
    }))
    fitP <- .nlmixr(.ocmt, d, "focei",
                    foceiControl(print = 0, maxOuterIterations = 2L,
                                 outerOpt = "nlminb", covMethod = "", calcTables = FALSE))

    # same columns as plain focei on the same model, natural theta order
    expect_identical(names(fitM$parHist), names(fitP$parHist))

    ph <- fitM$parHistData
    g <- ph[grepl("Gradient|Difference|Sensitivity", ph$type), ]
    expect_true(nrow(g) > 0)
    # mu columns: NaN gradients; optimizer columns: finite
    expect_true(all(is.nan(g$tka)))
    expect_true(all(is.nan(g$tcl)))
    expect_true(all(is.nan(g$tv)))
    expect_true(all(is.finite(g$add.sd)))
    expect_true(all(is.finite(g$o1)))

    s <- ph[ph$type == "Scaled", ]
    u <- ph[ph$type == "Unscaled", ]
    b <- ph[ph$type == "Back-Transformed", ]
    expect_true(nrow(u) > 0 && nrow(b) > 0)
    # mu columns pass through unscaled (# == U) and back-transform in X
    expect_identical(s$tka, u$tka)
    expect_equal(exp(u$tka), b$tka)
    expect_true(all(is.finite(u$tka)))
    # final recorded mu values are the fitted thetas
    expect_equal(unname(u$tcl[nrow(u)]), unname(fitM$theta[["tcl"]]),
                 tolerance = 1e-6)

    # console: standard columns, no bolt-on |mu| row, key note, blank grad cells
    expect_false(any(grepl("^\\|   mu\\|", out)))
    expect_true(any(grepl("mu-referenced thetas are regression-updated", out)))
    hdr <- grep("^\\|    #\\|", out, value = TRUE)[1]
    expect_false(is.na(hdr))
    .pos <- vapply(c("tka", "tcl", "tv", "add\\.sd"),
                   function(nm) as.numeric(regexpr(paste0("\\b", nm, "\\b"), hdr)),
                   numeric(1))
    expect_true(all(.pos > 0))
    expect_true(all(diff(.pos) > 0))
    gradRows <- grep("^\\|    [GFCMSA]\\|", out, value = TRUE)
    expect_true(length(gradRows) > 0)
    .blanks <- vapply(gradRows,
                      function(r) sum(gregexpr(" {11}\\|", r)[[1]] > 0),
                      numeric(1))
    expect_true(all(.blanks == 3))

    # plain focei is untouched: full-width columns, no NaN gradient entries
    phP <- fitP$parHistData
    expect_identical(ncol(phP), ncol(ph))
    gP <- phP[grepl("Gradient|Difference|Sensitivity", phP$type), ]
    expect_true(nrow(gP) > 0)
    expect_false(any(is.nan(gP$tka)))
  })

  test_that("mfocei parHist includes the mu covariate coefficient at its natural position", {
    skip_on_cran()

    theo_sd2 <- nlmixr2data::theo_sd
    theo_sd2$logWT <- log(theo_sd2$WT / 70)
    mod <- function() {
      ini({
        tka <- 0.45
        tcl <- 1
        tv <- 3.45
        allo.cl <- 0.75
        eta.ka ~ 0.6
        eta.cl ~ 0.3
        eta.v ~ 0.1
        add.sd <- 0.7
      })
      model({
        ka <- exp(tka + eta.ka)
        cl <- exp(tcl + eta.cl + allo.cl * logWT)
        v <- exp(tv + eta.v)
        linCmt() ~ add(add.sd)
      })
    }
    fit <- .nlmixr(mod, theo_sd2, "mfocei",
                   mfoceiControl(print = 0, maxOuterIterations = 2L,
                                 outerOpt = "nlminb", covMethod = "", calcTables = FALSE))
    nm <- names(fit$parHist)
    expect_true(all(c("tka", "tcl", "tv", "allo.cl", "add.sd") %in% nm))
    # natural ini order: coefficient between tv and add.sd
    expect_true(which(nm == "allo.cl") > which(nm == "tv"))
    expect_true(which(nm == "allo.cl") < which(nm == "add.sd"))
    g <- fit$parHistData[grepl("Gradient|Difference|Sensitivity",
                               fit$parHistData$type), ]
    expect_true(nrow(g) > 0)
    expect_true(all(is.nan(g$allo.cl)))
    expect_true(all(is.nan(g$tcl)))
    expect_true(all(is.finite(g$add.sd)))
  })
})
