#
#
# Data appear to be nearly 10 days of 1Hz data, which gives nearly 3/4 million
# points.  This can take a while to load.
#
#http://stackoverflow.com/questions/4756989/how-to-load-data-quickly-into-r
# subsequently:
# http://cran.r-project.org/web/packages/saves/index.html
#
fi <- "tmp.txt"
fi <- "B084.txt.gz"
system.time(dat <- read.table(fi, nrows=25e3, skip=17e3, header=TRUE))
#colClasses=c()
dat[dat==999999] <- NA
summary(dat)

#library(sqldf)
#system.time(dat3<-read.csv.sql(fi, sep="\t", nrows=86400, sql=sqlc))
