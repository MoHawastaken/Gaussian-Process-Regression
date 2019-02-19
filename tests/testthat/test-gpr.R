context("test-gpr")

#simple examples validated by calculation by hand

test_that("predict works", {
  X1 <-  matrix(c(-1/2, 1/2), nrow = 1)
  y1 <- c(4, 4)
  GPRobj1 <- GPR.polynomial$new(X1, y1, 1/2, 1/4, 1)
  expect_equivalent(GPRobj1$predict(0), c(2,1/8))
  
  # fstar = 1/3*(y_1 + y_2), Covstar = 1/3 for any X,y,xstar
  X2 <- matrix(c(1,2),nrow = 1)
  y2 <- c(1,3)
  GPRobj2 <- GPR.constant$new(X2, y2, 1, 1)
  expect_equivalent(GPRobj2$predict(3), c(4/3, 1/3))
  X3 <- matrix(c(100,54),nrow = 1)
  y3 <- c(5,0)
  GPRobj3 <- GPR.constant$new(X3, y3, 1, 1)
  expect_equivalent(GPRobj3$predict(pi), c(5/3, 1/3))

  # fstar = 1/(4-exp(-1)) c(2exp(-1/2)-exp(-5/2), 2exp(-2)-exp(-1)) %*% y
  #Covstar = cov for any y, with this choice of X,xstar
  X4 <- matrix(c(1, 2),nrow = 1)
  y4 = c(0, 1)
  GPRobj4 <- GPR.sqrexp$new(X4, y4, 1, 1)
  cov <- 1 - (2 * exp(-1) - 2 * exp(-3) + 2 * exp(-4)) / (4 - exp(-1))
  expect_equivalent(GPRobj4$predict(0), c((2*exp(-2) - exp(-1))/(4 - exp(-1)), cov))
})
