test_that("normalize works", {
  
  set.seed(1234)
  x <- rnorm(2000)
  
  pd <- spec.pgram(x, plot = FALSE, spans = 5)
  
  # spec test
  pd1a <- normalize(pd, 1, "spectrum", verbose = FALSE)
  pd1b <- normalize(pd, -1, "spectrum", verbose = FALSE)
  
  expect_error(normalize(pd, Fsamp = 0, src = 'spectrum', verbose = FALSE))
  expect_equal(pd1a, pd1b)

  # list test
  class(pd) <- NULL
  pd2 <- normalize(pd, Fsamp = 1, src = "spectrum", verbose = FALSE)
  class(pd2) <- 'spec'
  expect_equal(pd1a, pd2)
  
  
})
