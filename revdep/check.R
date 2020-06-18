library(devtools)
library(remotes)

if (!require(revdepcheck)){
    install_github('r-lib/revdepcheck')
    library(revdepcheck)
}

pkg <- ".."
revdepcheck::revdep_check(pkg, num_workers = 3)
revdep_report(pkg)
revdep_report_summary(pkg)
