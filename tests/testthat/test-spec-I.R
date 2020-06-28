
##

context("Spectrum estimation tools -- I")

tol <- 0.07

test_that("classes are correct",{
  
  set.seed(1234)
  x <- rnorm(100)
  pd <- stats::spectrum(x, plot=FALSE)
  pc <- psdcore(x, plot = FALSE, verbose = FALSE)
  pa <- pspectrum(x, plot = FALSE, verbose = FALSE)
  pa_b <- pspectrum_basic(x, verbose = FALSE)
  
  expect_is(pd, 'spec')
  expect_is(pc, c('amt','spec'))
  expect_is(pa, c('amt','spec'))
  expect_is(pa_b, c('amt','spec'))
  
  expect_is(normalize(pd, verbose = FALSE), 'spec')
  expect_is(normalize(pa, verbose = FALSE), c('amt','spec'))
  
  expect_is(pd[['taper']], 'numeric')
  expect_is(pa[['taper']], 'tapers')
  
})

test_that("pgram.compare results",{
  set.seed(1234)
  x <- rnorm(100)
  xp <- pspectrum(x, plot = FALSE, verbose=FALSE)
  expect_is(xp, 'amt')
  xpc <- pgram_compare(xp)
  expect_is(xpc, 'list')
})


test_that("pspectrum results are accurate",{
  
  set.seed(1234)
  x <- rnorm(100)
  varx <- var(x)
  twovar <- 2*varx
  
  xt <- ts(x, frequency=1)
  xt2 <- ts(x, frequency=10)
  
  pc <- pspectrum(xt, plot = FALSE, verbose = FALSE)
  pc2 <- pspectrum(xt2, plot = FALSE, verbose = FALSE)
  
  # make sure Nyquist frequencies are correct
  fn <- max(pc[['freq']])
  fn2 <- max(pc2[['freq']])
  expect_equal(fn, frequency(xt)/2)
  expect_equal(fn2, frequency(xt2)/2)
  
  # and normalization is for a single-sided spectrum
  nf <- pc[['numfreq']]
  nf2 <- pc2[['numfreq']]
  psum <- sum(pc[['spec']])
  psum2 <- sum(pc2[['spec']])
  expect_equal(fn*psum/nf, varx, tolerance=tol)
  expect_equal(fn2*psum2/nf2, varx, tolerance=tol)
  
  # normalization effects
  pcn <- normalize(pc, verbose = FALSE)
  pcn2 <- normalize(pc2, verbose = FALSE)
  nnf <- pcn[['numfreq']]
  nnf2 <- pcn2[['numfreq']]
  psumn <- sum(pcn[['spec']])
  psumn2 <- sum(pcn2[['spec']])
  expect_equal(fn*psumn/nnf, varx, tolerance=tol)
  expect_equal(fn2*psumn2/nnf2, varx, tolerance=tol)
  
})

test_that("psdcore arguments are tested",{
  set.seed(1234)
  x <- rnorm(100)
  xp1 <- psdcore(X.d = x, X.frq = 1, plot = FALSE)
  xp2 <- psdcore(X.d = x, X.frq = -1, plot = FALSE)
  expect_is(xp1, 'spec')
  expect_is(xp2, 'spec')
  expect_equal(xp1,xp2)
  expect_error(psdcore(X.d = x, X.frq = "1", plot = FALSE))
})

test_that("psdcore results are accurate",{
  
  set.seed(1234)
  x <- rnorm(100)
  varx <- var(x)
  twovar <- 2*varx
  
  xt <- ts(x, frequency=1)
  xt2 <- ts(x, frequency=10)
  
  pc <- psdcore(xt, verbose = FALSE, plot = FALSE)
  pc2 <- psdcore(xt2, verbose = FALSE, plot = FALSE)
  
  # make sure Nyquist frequencies are correct
  fn <- max(pc[['freq']])
  fn2 <- max(pc2[['freq']])
  expect_equal(fn, frequency(xt)/2)
  expect_equal(fn2, frequency(xt2)/2)
  
  # and normalization is for a single-sided spectrum
  nf <- pc[['numfreq']]
  nf2 <- pc2[['numfreq']]
  psum <- sum(pc[['spec']])
  psum2 <- sum(pc2[['spec']])
  expect_equal(psum/nf, twovar, tolerance=tol)
  expect_equal(psum2/nf2, twovar, tolerance=tol)
  
})

test_that("prewhiten results are accurate",{
  
  set.seed(1234)
  x <- rnorm(100)
  varx <- var(x)
  twovar <- 2*varx
  
  xt <- ts(x, frequency=1)
  xt2 <- ts(x, frequency=10)
  
  pc <- prewhiten(xt, plot = FALSE, verbose = FALSE)
  pc2 <- prewhiten(xt2, plot = FALSE, verbose = FALSE)
  
  # make sure Nyquist frequencies are correct
  fn <- frequency(pc[['prew_lm']])
  fn2 <- frequency(pc2[['prew_lm']])
  expect_equal(fn, frequency(xt))
  expect_equal(fn2, frequency(xt2))
  
  pa <- prewhiten(xt, AR.max = 10, plot = FALSE, verbose = FALSE)
  pa2 <- prewhiten(xt2, AR.max = 10, plot = FALSE, verbose = FALSE)
  
  # make sure Nyquist frequencies are correct
  fn <- frequency(pa[['prew_ar']])
  fn2 <- frequency(pa2[['prew_ar']])
  expect_equal(fn, frequency(xt))
  expect_equal(fn2, frequency(xt2))
  
})

test_that("pilot_spec results are accurate",{
  
  set.seed(1234)
  x <- rnorm(100)
  varx <- var(x)
  twovar <- 2*varx
  
  xt <- ts(x, frequency=1)
  xt2 <- ts(x, frequency=10)
  
  pc <- pilot_spec(xt, plot = FALSE, verbose = FALSE)
  pc2 <- pilot_spec(xt2, plot = FALSE, verbose = FALSE)
  
  # make sure Nyquist frequencies are correct
  fn <- max(pc[['freq']])
  fn2 <- max(pc2[['freq']])

  expect_equal(fn, frequency(xt)/2)
  expect_equal(fn2, frequency(xt2)/2)
  
  # and normalization is for a single-sided spectrum
  nf <- pc[['numfreq']]
  nf2 <- pc2[['numfreq']]
  psum <- sum(pc[['spec']])
  psum2 <- sum(pc2[['spec']])
  expect_equal(psum/nf, twovar, tolerance=tol)
  expect_equal(psum2/nf2, twovar, tolerance=tol)
  
  expect_warning(pa <- pilot_spec(xt, remove.AR = TRUE, plot = FALSE, verbose = FALSE)) # because there is no AR structure!
  expect_warning(pa2 <- pilot_spec(xt2, remove.AR = TRUE, plot = FALSE, verbose = FALSE))
  
  # make sure Nyquist frequencies are correct
  fn <- max(pa[['freq']])
  fn2 <- max(pa2[['freq']])
  
  expect_equal(fn, frequency(xt)/2)
  expect_equal(fn2, frequency(xt2)/2)
  
})


test_that("sampling rates for ts objects are honored",{
  set.seed(1234)
  x <- rnorm(100)
  f <- 1.2313
  xt <- stats::ts(x, frequency = f)
  
  p <- psdcore(xt)
  expect_equal(p[['nyquist.frequency']], f/2)
  expect_equal(f/2, max(p[['freq']]))
  
  pil <- pilot_spec(xt)
  expect_equal(pil[['nyquist.frequency']], f/2)
  expect_equal(f/2, max(pil[['freq']]))
  
  ps <- pspectrum(xt)
  expect_equal(ps[['nyquist.frequency']], f/2)
  expect_equal(f/2, max(ps[['freq']]))
  
  # make sure we cant trick ourselves into a new sampling rate
  expect_equal(p, psdcore(xt, X.frq = 99))
  expect_equal(pil, pilot_spec(xt, x.frequency = 99))
  expect_equal(ps, pspectrum(xt, x.frqsamp = 99))
  
  XT <- cbind(xt, xt)
  expect_equal(psdcore(XT), psdcore(XT, X.frq = 99))
  expect_equal(pilot_spec(XT), pilot_spec(XT, x.frequency = 99))
  expect_equal(pspectrum(XT), pspectrum(XT, x.frqsamp = 99))
})


test_that("check fast version",{
  set.seed(1234)
  x <- rnorm(100)
  expect_equal(psdcore(x, verbose = FALSE, plot = FALSE, fast = FALSE),
               psdcore(x, verbose = FALSE, plot = FALSE, fast = TRUE))
  #expect_equal(pspectrum(x, verbose = FALSE, plot = FALSE, fast = FALSE),
  #             pspectrum(x, verbose = FALSE, plot = FALSE, fast = TRUE))
  expect_equal(pilot_spec(x, verbose = FALSE, plot = FALSE, fast = FALSE),
               pilot_spec(x, verbose = FALSE, plot = FALSE, fast = TRUE))
  
  xt2 <- ts(x, frequency=10)
  expect_equal(psdcore(xt2, verbose = FALSE, plot = FALSE, fast = FALSE),
               psdcore(xt2, verbose = FALSE, plot = FALSE, fast = TRUE))
  #expect_equal(pspectrum(xt2, verbose = FALSE, plot = FALSE, fast = FALSE),
  #             pspectrum(xt2, verbose = FALSE, plot = FALSE, fast = TRUE))
  expect_equal(pilot_spec(xt2, verbose = FALSE, plot = FALSE, fast = FALSE),
               pilot_spec(xt2, verbose = FALSE, plot = FALSE, fast = TRUE))
})

test_that("check multivariate autospectra for psdcore",{

  set.seed(1234)
  x1 <- rnorm(100)
  x2 <- x1 * 2
  x <- cbind(x1, x2)  
  ps  <- psdcore(x, verbose = FALSE, plot = FALSE, fast = TRUE)
  
  ps1 <- psdcore(x[,1], verbose = FALSE, plot = FALSE, fast = TRUE)
  ps2 <- psdcore(x[,2], verbose = FALSE, plot = FALSE, fast = TRUE)
  
  expect_equal(ps$spec[,1], ps1$spec)
  expect_equal(ps$spec[,2], ps2$spec)
  
  
  ps1 <- psdcore(x[,1], verbose = FALSE, plot = FALSE, fast = FALSE)
  ps2 <- psdcore(x[,2], verbose = FALSE, plot = FALSE, fast = FALSE)
  
  expect_equal(ps$spec[,1], ps1$spec)
  expect_equal(ps$spec[,2], ps2$spec)
  
  # ps <- pspectrum(x)
})



test_that("check multivariate autospectra for pilot_spec",{
  
  set.seed(1234)
  x <- matrix(rnorm(200), ncol = 2)
  
  ps1 <- pilot_spec(x[,1], verbose = FALSE, plot = FALSE, fast = FALSE)
  ps2 <- pilot_spec(x[,2], verbose = FALSE, plot = FALSE, fast = FALSE)
  
  ps  <- pilot_spec(x, verbose = FALSE, plot = FALSE, fast = TRUE)
  
  expect_equal(ps$spec[,1], ps1$spec)
  expect_equal(ps$spec[,2], ps2$spec)
  
  psps <- pspectrum(x, plot = FALSE)
  
  expect_false(ps1[['is.multivariate']])
  expect_false(ps2[['is.multivariate']])
  expect_true(ps[['is.multivariate']])
  expect_true(psps[['is.multivariate']])
  
})


test_that("check multivariate autospectra, coherence, and phase for psdcore",{
  
  # compare spec.pgram to psdcore - the methods aren't equivalent but give 
  # similar results
  
  set.seed(1234)
  x <- matrix(rnorm(2000), ncol = 2)

  pd <- spec.pgram(x, plot = FALSE, spans = 5)
  pd <- normalize(pd, 1, "spectrum", verbose = FALSE)
  pc <- psdcore(x, plot = FALSE, verbose = FALSE,  ntaper = as.tapers(5))

  
  expect_null(pd[['is.multivariate']])
  expect_true(pc[['is.multivariate']])
  
  
  # is there an indexing issue here, need to remove the first value in psdcore
  # to get alignment

  pc$coh <- pc$coh[-1]
  pc$phase <- pc$phase[-1]
  pc$spec <- pc$spec[-1, ]
  
  n <- length(pd$coh)
  pd$coh <- pd$coh[-n]
  pd$phase <- pd$phase[-n]
  pd$spec <- pd$spec[-n, ]
  
  
  # select narrow range to avoid wrap around issues for phase
  sub_range <- 1.2
  wh1 <- which(pd$phase < pi/sub_range & pd$phase > -pi/sub_range)
  wh2 <- which(pc$phase < pi/sub_range & pc$phase > -pi/sub_range)
  wh  <- intersect(wh1, wh2)
  
  # do regression between methods
  phase_coef   <- coefficients(lm(pd$phase[wh]~pc$phase[wh]))
  coh_coef     <- coefficients(lm(pd$coh~pc$coh))
  spec_coef_1  <- coefficients(lm(pd$spec[,1]~pc$spec[,1]))
  spec_coef_2  <- coefficients(lm(pd$spec[,2]~pc$spec[,2]))
  
  
  # check intercept
  expect_lte(abs(phase_coef[1]), 0.05)
  expect_lte(abs(coh_coef[1]), 0.05)
  expect_lte(abs(spec_coef_1[1]), 0.05)
  expect_lte(abs(spec_coef_2[1]), 0.05)
  
  
  # check slope 
  expect_lte(abs(phase_coef[2]-1), 0.05)
  expect_lte(abs(coh_coef[2]-1), 0.05)
  expect_lte(abs(spec_coef_1[2]-1), 0.05)
  expect_lte(abs(spec_coef_2[2]-1), 0.05)

  
  # plot(pd$spec[,1])
  # points(pc$spec[,1], type='l')
  # plot(pd$spec[,2])
  # points(pc$spec[,2], type='l')
  # plot(pd$coh)
  # points(pc$coh, type='l')
  # plot(pd$phase)
  # points(pc$phase, type='l')

  
})

test_that("check multivariate output column selection works",{
  
  library(psd)
  data(wipp30)
  wipp30 <- as.matrix(wipp30[, -1])

  mv1 <- pspectrum(wipp30, riedsid_column= 0L, verbose = FALSE, plot = FALSE, output_column = 1)
  mv2 <- pspectrum(wipp30[, c(2,3,1)], riedsid_column= 0L, verbose = FALSE, plot = FALSE, output_column = 3)
  mv3 <- pspectrum(wipp30[, c(2,1,3)], riedsid_column= 0L, verbose = FALSE, plot = FALSE, output_column = 2)

  expect_equal(mv1, mv2)
  expect_equal(mv1, mv3)

  expect_true(mv1[['is.multivariate']])
  expect_true(mv2[['is.multivariate']])
  expect_true(mv3[['is.multivariate']])
  
})