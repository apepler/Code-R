rm(list=ls())
setwd("~/Documents/ECLs")
read.csv('mldbcomp.csv')->data

library(fields)
library("R.matlab")
readMat('~/Documents/Data/Useful.mat')->Useful
mask<-t(Useful$mask)
mask[is.na(mask)]=0
levels=seq(0,2,0.2)
cols=gray(seq(1,0.1,-0.1))

lats=seq(-49.5,-10.5,1)
lons=seq(110.5,179.5,1)
#lats=seq(-49.5,-10.5,2)
#lons=seq(110.5,179.5,2)


HR<-array(0,dim=c(length(lons),length(lats),4))
for(i in 1:length(lons))
  for(j in 1:length(lats))
  {
    I=which(data[,7]>=floor(lons[i]) & data[,7]<ceiling(lons[i]) & data[,6]>=floor(lats[j]) & data[,6]<ceiling(lats[j]))
    
    HR[i,j,1]=sum(data[I,14])/length(I) ##NCEP
    HR[i,j,2]=sum(data[I,15])/length(I) ##NCEP
    HR[i,j,3]=sum(data[I,17])/length(I) ##NCEP
    HR[i,j,4]=sum(data[I,19])/length(I) ##NCEP    
  }

image.plot(lons,lats,HR[,,1],xlab="",ylab="",breaks=seq(0,1,0.1),col=cols,zlim=c(0,1),xlim=c(145,161),ylim=c(-41,-24))
contour(Useful$x,Useful$y,mask,add=T,drawlabels=F)
image.plot(lons,lats,HR[,,4],xlab="",ylab="",breaks=seq(0,1,0.1),col=cols,zlim=c(0,1),xlim=c(145,161),ylim=c(-41,-24))
contour(Useful$x,Useful$y,mask,add=T,drawlabels=F)
image.plot(lons,lats,HR[,,1]-HR[,,4],xlab="",ylab="",breaks=seq(-0.5,0.5,0.1),col=cols,zlim=c(-0.5,0.5),xlim=c(145,161),ylim=c(-41,-24))
contour(Useful$x,Useful$y,mask,add=T,drawlabels=F)
