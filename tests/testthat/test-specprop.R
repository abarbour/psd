context("Spectral properties")

tap <- 1:10

sp1 <- spectral_properties(tap)
sp2 <- spectral_properties(as.tapers(tap))

set.seed(1234)
x <- rnorm(100)
pd <- stats::spectrum(x, plot=FALSE)
pc <- psdcore(x, plot = FALSE, verbose = FALSE)
pa <- pspectrum(x, plot = FALSE, verbose = FALSE)

sp3 <- spectral_properties(pd)
sp4 <- spectral_properties(pc)
sp5 <- spectral_properties(pa)

test_that("inheritance",{
  expected <- 'data.frame'
  expect_is(sp1, expected)
  expect_is(sp2, expected)
  expect_is(sp3, expected)
  expect_is(sp4, expected)
  expect_is(sp5, expected)
})

test_that("equality",{
  expect_equal(sp1, sp2)
})
##