# Clear memory
rm(list=ls(all=TRUE))

ssr <- function(){
 source("stevensRfunctions.r")
  source("hks.R")
  # Libraries used - uncomment below and comment above if you need to
  # source("http://idsource.googlecode.com/svn/trunk/R/stevensRfunctions.R")
  # source("http://idsource.googlecode.com/svn/trunk/R/hks.R")
}
ssr()

# Set filenames
#fnCHPDat <- "S/data/monthly_data_chp.csv"
fnCHPGP <- "supp_dataset_S02.csv"
fnStudDat <-"supp_dataset_S01.csv"
#fnStudDat <- "S/data/severity_data_28Nov2016.csv"
fnPopDat <- "Weekly_Pop_98-11_v4.csv"

# Load the CHP and GP data
#aDatCHP <- hks.read.clean.chp(fnCHPDat)
aDatGP <- hks.read.clean.gpw(fnCHPGP)
allDat <- hks.read.clean.study(fnStudDat)

# Subset the data
qDatVac <- allDat[allDat$qDat & allDat$vac,]
qDatNoVac <- allDat[allDat$qDat & allDat$novac,]
allDatVac <- allDat[allDat$vac,]
allDatNoVac <- allDat[allDat$novac,]

# Make a few variables required in the global scope
popDat <- read.csv(fnPopDat)
vecAg1HK <- as.vector(
  c(3/5*popDat[1,2]+sum(popDat[1,3:4])+4/5*popDat[1,5],
    1/5*popDat[1,5]+sum(popDat[1,6:10]),
	sum(popDat[1,11:14]),
	sum(popDat[1,15:19])))
assumedTotalSize <- sum(vecAg1HK)
vecAg1HK <- vecAg1HK/sum(vecAg1HK)
vecAg1Study <- table(qDatNoVac$ag1)
