# as.spec
# function (x, spans = NULL, kernel = NULL, taper = 0.1, pad = 0, 
#           fast = TRUE, demean = FALSE, detrend = TRUE, plot = TRUE, 
#           na.action = na.fail, ...)
# spg.out <- list(freq = freq, 
#                 spec = spec, 
#                 coh = NULL, 
#                 phase = NULL, 
#                 kernel = NULL, 
#                 df = Inf, # degrees freedom
#                 bandwidth = bandwidth, 
#                 n.used = N, 
#                 orig.n = N0, 
#                 series = series, 
#                 snames = colnames(x), 
#                 method = ifelse(TRUE, "Adaptive multitaper spectral density", "orig"), 
#                 taper = taper, 
#                 pad = 0, 
#                 detrend = TRUE, 
#                 demean = FALSE)
# class(spg.out) <- "spec"