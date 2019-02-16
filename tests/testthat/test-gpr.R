context("test-gpr")

X <-  matrix(c(-1/2, 1/2), nrow = 1)
y <- c(4, 4)
GPRobj <- GPR.polynomial$new(X, y, 1/4, 1, 1/2)


test_that("predict works", {
  expect_equivalent(GPRobj$predict(0), list(c(2), c(1/8)))
})
