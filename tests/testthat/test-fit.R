context("test-fit")
X <- matrix(seq(0, 1.1, by = 0.1), nrow = 1)

Y1 <- 3 * as.vector(X)
Y2 <- rep(5, 12)
Y3 <- 3 * as.vector(X) ^ 2 - 2 * as.vector(X)
Y4 <- 5 * exp(-as.vector(X) ^ 2)
Y5 <- 5 * exp(-as.vector(X) ^ 5)
Y6 <- 5 / (1 + as.vector(X) ^ 2)

test_that("fit", {
  expect_equivalent(fit(X,Y1,0.05,list("linear","constant","polynomial","sqrexp","gammaexp","rationalquadratic"))$cov,"linear")
  expect_equivalent(fit(X,Y2,0.05,list("linear","constant","polynomial","sqrexp","gammaexp","rationalquadratic"))$cov,"constant")
  expect_equivalent(fit(X,Y3,0.05,list("linear","constant","polynomial","sqrexp","gammaexp","rationalquadratic"))$cov,"polynomial")
  expect_equivalent(fit(X,Y4,0.05,list("linear","constant","polynomial","sqrexp","gammaexp","rationalquadratic"))$cov,"sqrexp")
  expect_equivalent(fit(X,Y5,0.05,list("linear","constant","polynomial","sqrexp","gammaexp","rationalquadratic"))$cov,"gammaexp")
  expect_equivalent(fit(X,Y6,0.05,list("linear","constant","polynomial","sqrexp","gammaexp","rationalquadratic"))$cov,"rationalquadratic")
})
