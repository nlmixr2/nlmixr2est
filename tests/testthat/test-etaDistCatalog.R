# The family codes the C++ dispatch switches on are lotri's ROW NUMBERS.
# .etaDistFamilyCode() returns which(.tab$name == nm) and src/etaDistFam.h
# hardcodes the answers (#define RXETADIST_GAMMA 13).  Nothing checks the two
# agree at run time: insert or reorder a row in lotri's catalog and every family
# at or after it silently rebinds -- dgamma becomes RXETADIST_INVGAMMA and the
# M-step maximizes the wrong density, with no error and a plausible-looking fit.
#
# So pin the whole table.  Adding a COLUMN to lotriEtaDists() is safe and these
# tests stay green; adding or moving a ROW breaks a test here rather than a fit
# in the field.  If a family is genuinely being added, append it and update both
# this list and the #define -- in the same commit.

.edCat <- function() {
  .t <- try(lotri::lotriEtaDists(), silent = TRUE)
  if (inherits(.t, "try-error")) testthat::skip("lotri::lotriEtaDists() unavailable")
  .t
}

test_that("the family name -> code table is exactly what the C++ #defines assume", {
  .t <- .edCat()
  .expect <- c(dnorm = 1L, stdNormal = 2L, studentT = 3L, dcauchy = 4L,
               doubleExponential = 5L, dlogis = 6L, gumbel = 7L, dlnorm = 8L,
               dchisq = 9L, invChiSquare = 10L, scaledInvChiSquare = 11L,
               dexp = 12L, dgamma = 13L, invGamma = 14L, dweibull = 15L,
               frechet = 16L, rayleigh = 17L, pareto = 18L, paretoType2 = 19L,
               dbeta = 20L, betaProportion = 21L, dunif = 22L)
  # the catalog itself has not shifted
  expect_identical(.t$name, names(.expect))
  # and the code each family resolves to is its position in that list
  .got <- vapply(names(.expect),
                 function(.n) .etaDistFamilyCode(str2lang(paste0(.n, "()"))),
                 integer(1))
  expect_identical(.got, .expect)
})

test_that("an unknown family takes code 0 and falls back to the R path", {
  expect_identical(.etaDistFamilyCode(str2lang("dNotAFamily()")), 0L)
})

test_that("roles come back in the family's own argument order", {
  .edCat()
  expect_identical(.etaDistRoles(str2lang("dgamma()")),
                   c("shape", "rate"))
  expect_identical(.etaDistRoles(str2lang("dnorm()")),
                   c("location", "scale"))
  expect_identical(.etaDistRoles(str2lang("studentT()")),
                   c("df", "location", "scale"))
  # rate is deliberately NOT folded into scale -- a covariate coefficient has
  # the opposite sign on the two, so grouping them would silently flip a slope
  expect_identical(.etaDistRoles(str2lang("dexp()")), "rate")
  expect_identical(.etaDistRoles(str2lang("rayleigh()")), "scale")
})

test_that("every family's roles line up with its parNames and are unique", {
  .t <- .edCat()
  if (is.null(.t$roles)) testthat::skip("lotri too old to carry roles")
  for (.i in seq_len(nrow(.t))) {
    .nm <- .t$name[.i]
    .roles <- .etaDistRoles(str2lang(paste0(.nm, "()")))
    .pars <- strsplit(.t$parNames[.i], ",", fixed = TRUE)[[1]]
    .pars <- .pars[nzchar(.pars)]
    expect_identical(length(.roles), length(.pars),
                     info = paste0(.nm, ": one role per argument"))
    # uniqueness is what makes a role a usable GROUP KEY: dbeta has
    # shape1/shape2, not two shapes, so "the shape group" is unambiguous
    expect_identical(anyDuplicated(.roles), 0L,
                     info = paste0(.nm, ": roles unique within the family"))
  }
})

test_that("role groups map back to argument positions", {
  .edCat()
  .g <- .etaDistRoleGroups(str2lang("dgamma()"))
  expect_identical(names(.g), c("shape", "rate"))
  expect_identical(.g$shape, 1L)
  expect_identical(.g$rate, 2L)
  .g3 <- .etaDistRoleGroups(str2lang("studentT()"))
  expect_identical(sort(unlist(.g3, use.names = FALSE)), 1:3)
})

test_that("unknown roles DECLINE rather than grouping everything together", {
  # the contract callers rely on: character(0)/NULL means "no role information",
  # which must not be read as "one big group"
  expect_identical(.etaDistRoles(str2lang("dNotAFamily()")),
                   character(0))
  expect_null(.etaDistRoleGroups(str2lang("dNotAFamily()")))
})

test_that("support endpoints are the roles that refuse a covariate", {
  .edCat()
  expect_true(all(c("lower", "upper") %in% .etaDistRoleNoCovariate))
  # dunif is bounds-only, so BOTH of its roles refuse a covariate
  expect_true(all(.etaDistRoles(str2lang("dunif()")) %in%
                    .etaDistRoleNoCovariate))
  # pareto's minimum is a support endpoint; its shape is not
  .p <- .etaDistRoles(str2lang("pareto()"))
  expect_true("lower" %in% .p)
  expect_false("shape" %in% .etaDistRoleNoCovariate)
})

test_that("support is reported for picking the warm-start surrogate", {
  .edCat()
  expect_identical(.etaDistSupport(str2lang("dgamma()")), "positive")
  expect_identical(.etaDistSupport(str2lang("dnorm()")), "real")
  expect_identical(.etaDistSupport(str2lang("dbeta()")), "unit")
  expect_identical(.etaDistSupport(str2lang("dNotAFamily()")),
                   NA_character_)
})
