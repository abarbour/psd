
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

test_that("riedsid2 works as expected",{
  expected <- 'integer'
  expect_is(riedsid2(pd), expected)
  expect_is(riedsid2(pc), expected)
  expect_is(riedsid2(pa), expected)
  expect_is(riedsid2(pa_b), expected)
})

test_that('riedsid issues deprecation warning',{
  expect_warning(riedsid(pd))
})


test_that("riedsid2 is equal to  riedsidrcpp",{
  
  expect_equal(riedsid2(pd), riedsid2(pd, fast = TRUE))
  expect_equal(riedsid2(pc), riedsid2(pc, fast = TRUE))
  expect_equal(riedsid2(pa), riedsid2(pa, fast = TRUE))
  expect_equal(riedsid2(pa_b), riedsid2(pa_b, fast = TRUE))
  
})
