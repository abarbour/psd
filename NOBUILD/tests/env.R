rm(list=ls(all=TRUE))
ls(all=T)

e <- new.env()
ls(all=T)

FUN1 <- function(...) new.env( parent=parent.frame(1) )
ef <- FUN1()
ls(all=T)

FUN2 <- function(...) new.env( parent=parent.frame(2) )
efp <- FUN2()
ls(all=T)

FUNA1 <- function(...) assign("xe", FUN1())
FUNA1()
ls(all=T)

FUNA2 <- function(...) assign("xe", FUN1(), pos=1)
FUNA2()
ls(all=T)

FUNA3 <- function(...) assign("xe", FUN1(), pos="somepos")
# DOESNT WORK: FUNA3()

FUNG1 <- function(...) get("xe", pos=1)
FUNG1()
ls(all=T)

envName = "someenv"
as.symbol(envName) = new.env()
ls(all=T)
