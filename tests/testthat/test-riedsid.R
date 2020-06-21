library(testthat)
##

context("Optimal taper estimation")

set.seed(1234)
x <- rnorm(200)

pd <- stats::spectrum(x, plot=FALSE)
pc <- psdcore(x, plot = FALSE, verbose = FALSE)
pa <- pspectrum(x, plot = FALSE, verbose = FALSE)
pa_b <- pspectrum_basic(x, verbose = FALSE)

FIG <- function(){
  plot(normalize(pd))
  plot(pc, add=TRUE, lwd=2)
  plot(pa, add=TRUE, col='red', lwd=2)
  plot(pa_b, add=TRUE, col='red', lty=2)
}

test_that("riedsid2 returns integer as expected",{
  expected <- 'integer'
  expect_is(riedsid2(pd), expected)
  expect_is(riedsid2(pc), expected)
  expect_is(riedsid2(pa), expected)
  expect_is(riedsid2(pa_b), expected)
})

test_that('riedsid issues deprecation warning in favor of riedsid2',{
  expect_warning(riedsid(pd))
})

test_that("riedsid2 R-version is equal to Rcpp version",{
  
  expect_equal(riedsid2(pd, fast=FALSE), riedsid2(pd, fast = TRUE))
  expect_equal(riedsid2(pc, fast=FALSE), riedsid2(pc, fast = TRUE))
  expect_equal(riedsid2(pa, fast=FALSE), riedsid2(pa, fast = TRUE))
  expect_equal(riedsid2(pa_b, fast=FALSE), riedsid2(pa_b, fast = TRUE))
  
})


test_that("multivariate riedsid2 works",{
  
  set.seed(1234)
  x <- matrix(rnorm(200), ncol = 2)
  taps <- ceiling(runif(200/2, 10, 300))
  
  pd <- stats::spectrum(x, plot=FALSE)
  
  # each separately and then take the minimum number of tapers
  r_s <- cbind(riedsid2(pd$spec[, 1], fast=FALSE),
               riedsid2(pd$spec[, 2], fast=FALSE))
  r_s <- apply(r_s, 1, min)
  
  # multivariate method
  r_mv <- riedsid2(pd$spec, fast=TRUE)
  expect_equal(r_mv, r_s)
  
  # spec method works
  r_mv_spec <- riedsid2(pd, fast=TRUE)
  expect_equal(r_mv_spec, r_s)
  
  
})


test_that("riedsid_rcpp  work",{
  set.seed(1234)
  x <- matrix(rnorm(200), ncol = 2)
  pd <- stats::spectrum(x, plot=FALSE)
  
  r_s1<- riedsid_rcpp(PSD = as.matrix(pd$spec[,1]), ntaper = 3, riedsid_column = 0)
  r_s2<- riedsid_rcpp(PSD = as.matrix(pd$spec[,1]), ntaper = 3, riedsid_column = -1)
  r_s3<- riedsid_rcpp(PSD = as.matrix(pd$spec[,1]), ntaper = 3, riedsid_column = 1)

  expect_equal(r_s1, r_s2)
  expect_equal(r_s2, r_s3)
  expect_warning(riedsid_rcpp(PSD = as.matrix(pd$spec[,1]), 
                            ntaper = 3, 
                            riedsid_column = 2))
  
})
  
  