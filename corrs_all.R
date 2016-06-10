rm(list=ls())
setwd('/media/Seagate Expansion Drive/monthly rainfall')
source('~/Documents/R/load.AWAP.season.ff.R')
library(ff)
m1<-c("06","07","08","09","10")
m2<-seq(6,10,1)
name<-"JJASO"
rain<-getRainff(m1)

# rm(list=ls())
# setwd('G:/monthly-rainfall')
# source('~/Documents/R/load.AWAP.wseason.ff.R')
# library(ff)
# m1a<-c("11","12")
# m1b<-c("01","02","03")
# m2a<-seq(11,12,1)
# m2b<-seq(1,3,1)
# name<-"NDJFM"
# rain<-getRainWff(m1a,m1b)

read.table('~/Documents/Timeseries/dmi.txt',header=T,sep="")->dmi
DMI=rowMeans(dmi[dmi[,1]>=1900,(m2[1]+1):(m2[length(m2)]+1)])
read.table('~/Documents/Timeseries/n34.txt',header=T,sep="")->n34
N34=rowMeans(n34[n34[,1]>=1900,(m2[1]+1):(m2[length(m2)]+1)])
read.table('~/Documents/Timeseries/iodE.txt',header=T,sep="")->iodE
IODE=rowMeans(iodE[iodE[,1]>=1900,(m2[1]+1):(m2[length(m2)]+1)])
read.table('~/Documents/Timeseries/iodW.txt',header=T,sep="")->iodW
IODW=rowMeans(iodW[iodW[,1]>=1900,(m2[1]+1):(m2[length(m2)]+1)])

# GDI=rowMeans(cbind(gdi[(length(gdi[,1])-111):(length(gdi[,1])-1),(m2a[1]+1):(m2a[length(m2a)]+1)],gdi[(length(gdi[,1])-110):(length(gdi[,1])),(m2b[1]+1):(m2b[length(m2b)]+1)]))

corrE<-corrW<-corrE_N<-corrW_N<-matrix(0,691,886)
source('~/Documents/R/pcor.R')

for(i in 1:691)
  for(j in 1:886)
  {
    corrE[i,j]<-cor(IODE,rain[i,j,])
    corrW[i,j]<-cor(IODW,rain[i,j,])
    a<-pcor.test(IODE,rain[i,j,],N34)
    b<-as.numeric(a[1])
    corrE_N[i,j]<-b
    a<-pcor.test(IODW,rain[i,j,],N34)
    b<-as.numeric(a[1])
    corrW_N[i,j]<-b
  }

library("R.matlab")
readMat('~/Documents/Data/Useful.mat')->Useful
mask<-t(Useful$mask)
mask[is.na(mask)]=0

source('~/Documents/R/color.palette.R')
pal <- color.palette(c("red","yellow","white","cyan","blue"), c(20,5,5,20))
cm=pal(20)
cm[9:12]="white"

source('~/Documents/R/corr.plot.R')
corr.plot(Useful$x,Useful$y,t(corrW_N*Useful$mask),mask,cm)
corr.plot(Useful$x,Useful$y,t(corrE*Useful$mask),mask,cm,paste('~/Documents/corrIODE_',name,'.pdf',sep=""))
