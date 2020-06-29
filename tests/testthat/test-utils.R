
##

context("Utility functions")

test_that("na_locf expected behavior", {
  n <- 5
  ntot <- 2 * n + 1
  y0 <- y0na <- seq_len(ntot)
  y0na[n + 1] <- NA
  numna <- 3
  nas <- rep(NA, numna)
  yl <- c(nas, y0)
  yr <- c(y0, nas)
  ylr <- c(nas, y0, nas)
  ylrna <- c(nas, y0na, nas)
  bad <- NULL
  
  expect_is(y0, 'integer')
  expect_is(na_locf(y0), 'numeric')
  expect_equal(y0, na_locf(y0))
  expect_equal(y0na, na_locf(y0na))
  expect_equal(as.numeric(as.character(y0na)), na_locf(y0na))
  
  ylLocf <- na_locf(yl)
  yrLocf <- na_locf(yr)
  ylrLocf <- na_locf(ylr)
  ylrnaLocf <- na_locf(ylrna)
  expect_length(ylLocf, length(yl))
  expect_length(yrLocf, length(yr))
  expect_length(ylrLocf, length(ylr))
  expect_length(ylrnaLocf, length(ylrna))
  expect_equal(range(ylLocf, na.rm=TRUE), range(yl, na.rm=TRUE))
  expect_equal(range(yrLocf, na.rm=TRUE), range(yr, na.rm=TRUE))
  expect_equal(range(ylrLocf, na.rm=TRUE), range(ylr, na.rm=TRUE))
  expect_equal(range(ylrnaLocf, na.rm=TRUE), range(ylrna, na.rm=TRUE))
  
  expect_warning(expect_identical(bad, na_locf(bad)))
  
  Y <- cbind(yl, yr)
  expect_equal(dim(Y), dim(na_locf(Y)))
  expect_true(!any(is.na(na_locf(Y))))
})

test_that("message delivery is working",{
  expect_message(adapt_message(0))
  expect_is(adapt_message(0), 'character')
  expect_message(adapt_message(1))
  expect_error(adapt_message(-1))
})

test_that("decibel conversions are accurate",{
  
  expect_equal(dB(2), 10*log10(2))
  expect_equal(dB(1), 0)
  # from dB to amp or power
  expect_equal(dB(0, invert=TRUE), 1)
  expect_equal(dB(0, invert=TRUE, is.power = TRUE), 1)
  # power ratio
  expect_equal(dB(70, invert = TRUE, is.power = FALSE), 1e7)
  # amplitude ratio
  expect_equal(round(dB(70, invert = TRUE, is.power = TRUE)), 3162)
})

test_that("variance of difference series are accurate",{
  X <- 1:10
  expect_equal(vardiff(X), 0)
  expect_equal(varddiff(X), 0)
})

test_that("polygon creation is accurate",{
  nx <- 10
  X <- seq_len(nx)
  Y <- X^2
  YE <- runif(nx, 2, 5)
  
  PX <- create_poly(X, Y, YE, from.lower=FALSE)
  PXl <- create_poly(X, Y, YE, from.lower=TRUE)
  
  expect_is(PX, 'data.frame')
  expect_is(PXl, 'data.frame')
  
  expect_equal(unique(PX$x.x), X)
  expect_equal(unique(PXl$x.x), X)
  
  expect_identical(PX$x.x, PXl$x.x)

})

test_that("colvec reshapes correctly", {
  m <- na_mat(2,2)
  expect_is(colvec(m), 'matrix')
  expect_equal(dim(colvec(m)), c(4,1))
})

test_that("rowvec reshapes correctly", {
  m <- na_mat(2,2)
  expect_is(rowvec(m), 'matrix')
  expect_equal(dim(rowvec(m)), c(1,4))
})

test_that("NA matrices are assembled correctly", {
  
  expect_is(na_mat(1),'matrix')
  expect_error(na_mat(-1))
  expect_error(na_mat(1,-1))
  
  expect_equal(dim(na_mat(2)), c(2,1))
  expect_equal(dim(na_mat(2,2)), c(2,2))
  expect_equal(dim(na.omit(na_mat(2,2))), c(0,2))
  
  expect_equal(range(na_mat(2,2)), c(NA_real_, NA_real_))
  expect_warning(range(na_mat(2,2), na.rm=TRUE))
               
})
  
test_that("zeros and ones are assembled correctly", {
  
  expect_is(zeros(1),'matrix')
  expect_equal(dim(zeros(2)), c(2,1))
  expect_equal(range(zeros(2)), c(0,0))
  
  expect_is(ones(1),'matrix')
  expect_equal(dim(ones(2)), c(2,1))
  expect_equal(range(ones(2)), c(1,1))
  
})

test_that("modulo division values are accurate", {

  expect_equal(mod(1,2), 1%%2)
  expect_equal(mod(1+1,2+1), (1+1)%%(2+1))
  
  mx <- 3
  X <- -mx:mx
  
  expect_equal(modulo_floor(X), modulo_floor(X, 2))
  expect_equal(modulo_floor(X,-2), modulo_floor(X, 2))
  
  n <- 1; expect_equal(range(modulo_floor(X, n)), range(X))
  n <- mx; expect_equal(range(modulo_floor(X, n)), n*c(-1,1))
  n <- mx+1; expect_equal(range(modulo_floor(X, n)), c(0,0))
  expect_equal(range(modulo_floor(X, NA)), c(0,0))
  
})

test_that("error checking works", {
  
  expect_warning(modulo_floor(1:10, Inf))
  
  expect_error(modulo_floor(1:10,0))
  expect_error(modulo_floor(1:10,NULL))
  expect_error(modulo_floor(1:10,1:2))
  expect_error(modulo_floor(1:10,-2:-1))
  
})

test_that('plotting functions return correctly', {
  set.seed(1234)
  x <- rnorm(100)
  pc <- psdcore(x, plot = FALSE, verbose = FALSE)
  expect_equal(plot(pc), lines(pc))
})


test_that('vector_reshape works', {
  
  x <- matrix(rnorm(10), ncol = 1)
  
  expect_equal(vector_reshape(x, vec.shape = 'horizontal'),t(x))
  expect_equal(vector_reshape(x, vec.shape = 'vertical'),as.matrix(x))
  
})

test_that('checking classes works', {
  s <- pspectrum(rnorm(10))
  
  expect_equal(is.spec(s), TRUE)
  expect_equal(is.amt(s), TRUE)
  expect_equal(is.tapers(s$taper), TRUE)

})

##
