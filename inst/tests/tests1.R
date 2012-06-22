# testfile 1

load("compareWITH.RData")



test_that("output of 3pl mle and wle works fine", {
  
  expect_that(PP_3PLwle(u=c(1,0,1,0),a=c(1.1,2,0.2,1),s=c(-2,-1.23,1.11,3.48),i=c(0,0.1,0.3,0.2)), equals(erg1))
  expect_that(PP_3PLwle(u=c(0,0,0,0),a=c(1.1,2,0.2,1),s=c(-2,-1.23,1.11,3.48),i=c(0,0.1,0.3,0.2)), equals(erg2))
  
  expect_that(PP_3PLmle(u=c(1,0,1,0),a=c(1.1,2,0.2,1),s=c(-2,-1.23,1.11,3.48),i=c(0,0.1,0.3,0.2)), equals(erg3))
  expect_that(PP_3PLmle(u=c(0,0,0,0),a=c(1.1,2,0.2,1),s=c(-2,-1.23,1.11,3.48),i=c(0,0.1,0.3,0.2)), equals(erg4))
  
})


test_that("SUMMARY of 3pl/GPCM mle and wle works fine", {
  
  expect_that(summary(PP_3PLwle(u=c(1,0,1,0),a=c(1.1,2,0.2,1),s=c(-2,-1.23,1.11,3.48),i=c(0,0.1,0.3,0.2))), equals(summary(erg1)))
  expect_that(summary(PP_3PLwle(u=c(0,0,0,0),a=c(1.1,2,0.2,1),s=c(-2,-1.23,1.11,3.48),i=c(0,0.1,0.3,0.2))), equals(summary(erg2)))
  
  expect_that(summary(PP_3PLmle(u=c(1,0,1,0),a=c(1.1,2,0.2,1),s=c(-2,-1.23,1.11,3.48),i=c(0,0.1,0.3,0.2))), equals(summary(erg3)))
  expect_that(summary(PP_3PLmle(u=c(0,0,0,0),a=c(1.1,2,0.2,1),s=c(-2,-1.23,1.11,3.48),i=c(0,0.1,0.3,0.2))), equals(summary(erg4)))
  
})



test_that("all random of 3pl/GPCM mle and wle works fine", {
  set.seed(10101)
  u_1 <- sample(0:1,15,replace=TRUE)
  a_1 <- runif(15, 0.4, 3)
  s_1 <- rnorm(15, 0, 2)
  i_1 <- runif(15, 0, 0.4)
  
  expect_that(PP_3PLmle(u=u_1,a=a_1,s=s_1,i=i_1), equals(ergZF1))
  expect_that(PP_3PLmle(u=u_1,a=a_1,s=s_1,i=i_1), equals(ergZF2))
  
})



test_that("Rasch Model functions", {
  pam1 <- c(-4.4,-3.448,0.1,0.2,1.11,2.18,3.01) 
  
  expect_that(PP_RMmle(pam1,expol=TRUE)  , equals(ergrm1))
  expect_that(PP_RMmle(pam1,expol=FALSE)  , equals(ergrm2))
  expect_that(PP_RMwle(pam1), equals(ergrm3))
  
  
  expect_that(summary(PP_RMmle(pam1,expol=TRUE))  , equals(summary(ergrm1)))
  expect_that(summary(PP_RMmle(pam1,expol=FALSE))  , equals(summary(ergrm2)))
  expect_that(summary(PP_RMwle(pam1)), equals(summary(ergrm3)))
})



test_that("GPCM (wle / mle) works fine", {
  
  THRES  <- matrix(c(-2.41,0.36,-0.97,0.92,-0.65,-0.65,-0.28,0.72),nrow=2)
  sl     <- c(0.5,1,1.5,1.1)
  v1     <- c(1,2,1,3)
  
  THRES2  <- matrix(c(-2.41,NA,-0.57,0.92,-0.65,-0.65,-0.28,1.19),nrow=2)
  sl2     <- c(1,1,1,1)
  v2     <- c(1,2,1,3)

  expect_that(PP_GPCMmle(u=v1,sl=sl,s=THRES), equals(erg5))
  expect_that(PP_GPCMwle(u=v1,sl=sl,s=THRES), equals(erg6))
  
  expect_that(PP_GPCMmle(u=v2,sl=sl2,s=THRES2), equals(erg7))
  expect_that(PP_GPCMwle(u=v2,sl=sl2,s=THRES2), equals(erg8))
})



test_that("SUMMARY of GPCM (wle / mle) works fine", {
  
  THRES  <- matrix(c(-2.41,0.36,-0.97,0.92,-0.65,-0.65,-0.28,0.72),nrow=2)
  sl     <- c(0.5,1,1.5,1.1)
  v1     <- c(1,2,1,3)
  
  THRES2  <- matrix(c(-2.41,NA,-0.57,0.92,-0.65,-0.65,-0.28,1.19),nrow=2)
  sl2     <- c(1,1,1,1)
  v2     <- c(1,2,1,3)
  
  expect_that(summary(PP_GPCMmle(u=v1,sl=sl,s=THRES)), equals(summary(erg5)))
  expect_that(summary(PP_GPCMwle(u=v1,sl=sl,s=THRES)), equals(summary(erg6)))
  
  expect_that(summary(PP_GPCMmle(u=v2,sl=sl2,s=THRES2)), equals(summary(erg7)))
  expect_that(summary(PP_GPCMwle(u=v2,sl=sl2,s=THRES2)), equals(summary(erg8)))
})


test_that("errors 3PL", {
  
  expect_that(PP_3PLwle(u=c(1,0,1,0,1),a=c(1.1,2,0.2,1),s=c(-2,-1.23,1.11,3.48),i=c(0,0.1,0.3,0.2)), throws_error())
  expect_that(PP_3PLwle(u=c(2,0,1,0),a=c(1.1,2,0.2,1),s=c(-2,-1.23,1.11,3.48),i=c(0,0.1,0.3,0.2)), throws_error())
  expect_that(PP_3PLwle(u=c(1,1,1,0),a=c(1.1,2,0.2,1),s=c(-2,-1.23,1.11,3.48),i=c(0,0.1,0.3,0.2,1)), throws_error())
  expect_that(PP_3PLwle(u=c(1,1,1,0),a=c(1.1,2,0.2,1),s=c(-2,-1.23,1.11,3.48,1),i=c(0,0.1,0.3,0.2)), throws_error())
})



test_that("errors GPCM", {
  
  THRES  <- matrix(c(-2.41,0.36,-0.97,0.92,-0.65,-0.65,-0.28,0.72),nrow=2)
  sl     <- c(0.5,1,1.5,1.1)
  v1     <- c(1,2,1,3)
  
  THRES2  <- matrix(c(-2.41,NA,-0.57,0.92,-0.65,-0.65,-0.28,1.19),nrow=2)
  sl2     <- c(1,1,1,1)
  v2     <- c(1,2,1,3)
  
  expect_that(PP_GPCMmle(u=c(1,2,1,3,1),sl=sl,s=THRES), throws_error()) # error when u = too long
  expect_that(PP_GPCMmle(u=c(1,2),sl=sl,s=THRES), throws_error()) # error when u = too short
  
  expect_that(PP_GPCMmle(u=v1,sl=c(0.5,1,1.5,1.1,1),s=THRES), throws_error()) # error when sl = too long
  expect_that(PP_GPCMmle(u=v1,sl=c(0.5),s=THRES), throws_error()) # error when sl = too short
  
  
})





