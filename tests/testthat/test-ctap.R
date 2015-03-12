
context("Taper constraints -- minspan")

test_that("Length and positivity requirements",{
  expect_error(minspan(1))
  expect_error(minspan(0))
  expect_error(minspan(-1))
  expect_error(minspan(-1:0))
})

test_that("Number of tapers limited by section length", {
  
  ms. <- minspan(0:2)
  expect_equal(min(ms.), 1)
  expect_equal(max(ms.), 1)
  
  ms. <- minspan(0:3)
  expect_equal(min(ms.), 1)
  expect_equal(max(ms.), 2)
  
  ms. <- minspan(0:4)
  expect_equal(min(ms.), 1)
  expect_equal(max(ms.), 2)
  
  ms. <- minspan(0:5)
  expect_equal(min(ms.), 1)
  expect_equal(max(ms.), 3)
  
  ms. <- minspan(0:6)
  expect_equal(min(ms.), 1)
  expect_equal(max(ms.), 3)
  
  ms. <- minspan(0:7)
  expect_equal(min(ms.), 1)
  expect_equal(max(ms.), 4)
})

test_that("Strange values are dealt with", {
  
  expect_warning(ms. <- minspan(c(0:7,Inf)))
  expect_equal(min(ms.), 1)
  expect_equal(max(ms.), 4)
  expect_equal(length(ms.), 9)
  
  ms. <- minspan(c(0:7,NA))
  expect_equal(min(ms.), 1)
  expect_equal(max(ms.), 4)
  expect_equal(length(ms.), 9)
  
  ms. <- minspan(c(0:7,""))
  expect_equal(min(ms.), 1)
  expect_equal(max(ms.), 4)
  expect_equal(length(ms.), 9)
  
  ms. <- minspan(c(0:7,NULL))
  expect_equal(min(ms.), 1)
  expect_equal(max(ms.), 4)
  expect_equal(length(ms.), 8) # instead of 9
  
})
