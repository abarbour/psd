
if (!require(revdepcheck)){
    library(remotes)
    install_github('r-lib/revdepcheck')
    library(revdepcheck)
}

pkg <- ".."
revdepcheck::revdep_check(pkg, num_workers = 3)
revdepcheck::revdep_report(pkg)
revdepcheck::revdep_report_summary(pkg)
