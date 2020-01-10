test_that("coherence and phase works for single entry", {
  
  a <- array(rnorm(100), dim = c(100, 1, 1))
  
  expect_warning(phase(a))
  expect_warning(coherence(a))
  
  
})
