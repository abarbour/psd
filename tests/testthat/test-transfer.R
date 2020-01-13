test_that("transfer function works", {
  
  # can we apply and recover a constant response
  set.seed(1234)
  seq_length  <- 10000
  x <- sin(seq(0, 8*pi, length.out = seq_length)) + rnorm(seq_length)
  y <- x * 0.5
  xy <- cbind(y, x)
  pc <- psdcore(xy, plot = FALSE, verbose = FALSE,  ntaper = as.tapers(3))
  
  expect_equal(Mod(pc$transfer), 
               matrix(0.5, nrow = length(pc$transfer), ncol = 1))
  
  
  # can we apply and recover a variable response
  # generate data
  set.seed(1234)
  kern_length <- 101
  seq_length  <- 10000
  x <- sin(seq(0, 8*pi, length.out = seq_length)) + rnorm(seq_length)
  kern <- exp(seq(-10, 0, length.out = kern_length)) / kern_length
  pad_kern <- c(rep(0, seq_length-kern_length), kern)
  # convolution
  y <- Re(fft(fft(x) * Conj(fft(pad_kern)), inverse = TRUE)) / seq_length
  
  # plot(x, type='l')
  # points(y, col = 'red', type='l')

  xy <- cbind(y, x)
  pc <- psdcore(xy, plot = FALSE, verbose = FALSE,  ntaper = as.tapers(3))
  
  # convert to impulse response
  # this is necessary for the inverse fft as we use half length
  pc$transfer <- c(pc$transfer, rev(Conj(pc$transfer)))
  
  # inverse fft of transfer fun
  n <- floor(NROW(pc$transfer) / 2) + 1
  imp <- fft(pc$transfer, inverse = TRUE) / NROW(pc$transfer)
  imp <- Mod(imp) * sign(Re(imp))
  mt_impulse <- cumsum(rev(imp)[1:n] + (imp)[1:n])
  
  expect_equal(cumsum(rev(kern)), mt_impulse[1:kern_length], tolerance = 1e-3)
  
  # plot(cumsum(rev(kern)), type='l', col = 'red')
  # points(mt_impulse[1:kern_length], cex = 0.5, pch = 20)
  
  
})
