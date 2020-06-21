test_that("parabolic weight methods give equivalent results", {
  
  w1 <- parabolic_weights_rcpp(ntap = 5)$taper_weights
  w2 <- as.numeric(parabolic_weights_field(ntap = 5)[[5]])
  
  expect_equal(w1, w2)

})
