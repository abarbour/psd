if (!require(revdepcheck)) source("https://install-github.me/r-lib/revdepcheck") #install
library(revdepcheck)
res <- revdepcheck::revdep_check(num_workers = 4)
revdepcheck::revdep_check()

library(withr)
withr::with_output_sink(
  "revdep/cran.md",
  revdepcheck::revdep_report_cran()
)
