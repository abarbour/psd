
##

context("Spectrum estimation tools -- II")

test_that("'updating' spectrum is not yet implemented",{
  set.seed(1234)
  x <- rnorm(100)
  
  pd <- stats::spectrum(x, plot=FALSE)
  expect_error(pspectrum(pd))
  
  pc <- psdcore(x, plot = FALSE, verbose = FALSE)
  expect_error(pspectrum(pc))
  
  pa <- pspectrum(x, plot = FALSE, verbose = FALSE)
  expect_error(pspectrum(pa))
  
  pa_b <- pspectrum_basic(x, verbose = FALSE)
  expect_error(pspectrum(pa_b))
})
