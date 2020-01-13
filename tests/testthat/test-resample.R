test_that("resampling methods give same results", {
  
  n. <- 100000
  set.seed(1234)
  nc <- 2
  x <- matrix(cumsum(sample(c(-1, 1), n., TRUE)), ncol=nc)
  fftz <- mvfft(x)
  psd <- Re(fftz * Conj(fftz))
  taps <- ceiling(runif(n./nc,10,300))
  
  rsz1a <- resample_fft_rcpp(fftz[,1], taps, verbose = FALSE)
  rsz1b <- resample_fft_rcpp2(fftz[,1], taps, verbose = FALSE)
  
  rsz2a <- resample_fft_rcpp(fftz[,2], taps, verbose = FALSE)
  rsz2b <- resample_fft_rcpp2(fftz[,2], taps, verbose = FALSE)
  
  rsz3 <- resample_mvfft(fftz, taps, verbose = FALSE)
  

  expect_equal(as.numeric(rsz1a$psd), as.numeric(Re(rsz3$psd[,1,1])))
  expect_equal(as.numeric(rsz2a$psd), as.numeric(Re(rsz3$psd[,2,2])))
  
  expect_equal(as.numeric(rsz1b$psd), as.numeric(Re(rsz3$psd[,1,1])))
  expect_equal(as.numeric(rsz2b$psd), as.numeric(Re(rsz3$psd[,2,2])))
  
})

test_that("riedsid_rcpp gives minimum of matrix", {
  
  n. <- 100000
  set.seed(1234)
  nc <- 2
  x <- matrix(cumsum(sample(c(-1, 1), n., TRUE)), ncol=nc)
  fftz <- mvfft(x)
  psd <- Re(fftz * Conj(fftz))
  taps <- ceiling(runif(n./nc,10,300))
  
  
  a1 <- riedsid_rcpp(as.matrix(psd)[, 1, drop = FALSE], taps)
  a2 <- riedsid_rcpp(psd[, 2, drop = FALSE], taps)
  b  <- riedsid_rcpp(as.matrix(psd), taps)
  
  expect_equal(apply(cbind(a1,a2), 1, min), as.numeric(b))
  
})
