
context("Discrete Fourier Transform calculations")

library(stats)
stopifnot(require(fftw))

n. <- 10
x. <- seq_len(n.)
xn. <- 1.0*x.

test_that("fftw expects numeric or complex", {
  expect_error(FFT(x.))
})

test_that("fft and fftw return complex", {
  expect_is(fft(xn.), 'complex')
  expect_is(FFT(xn.), 'complex')
})

test_that("fft and fftw return equivalent results", {
  
  # Forward transform
  expect_equal(fft(x.), FFT(xn.))
  expect_equal(fft(xn.), FFT(xn.))
  
  # Inverse transform
  # - by default FFTW scales the inverse transform so this is an error:
  expect_error(stopifnot(all.equal(fft(xn., inverse = TRUE), IFFT(xn.))))
  # but this is not:
  expect_equal(fft(xn., inverse = TRUE), IFFT(xn., scale=FALSE))
  
})
