#' \dontrun{
##
## modulo division
##
x <- 1:10
mc1a <- mod(1,2)
mc2a <- mod(1+x,2)
mc1b <- 1 %% 2
mc2b <- 1 + x %% 2
mc2c <- (1 + x) %% 2
all.equal(mc1a, mc1b) # TRUE
all.equal(mc2a, mc2b) # "Mean absolute difference: 2"
all.equal(mc2a, mc2c) # TRUE
##
#' }
