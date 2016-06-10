rm(list=ls())
setwd('~/Documents/ECLs')
ncep1<-read.csv('CSV/ECLfixes_NCEP1.csv')
lats=seq(-40.5,-23.5,1)
lons=seq(149.5,160.5,1)
library(fields)
library("R.matlab")
readMat('~/Documents/Data/Useful.mat')->Useful
mask<-t(Useful$mask)
mask[is.na(mask)]=0
years=seq(1950,2012)
rain<-read.csv('morwenna.csv')
rain2=rain[which(rain[,6]==1),] 
ncep1[,8]=ncep1[,1] %in% rain2[,3]
yy=floor(ncep1[,1]/10000)
cols=gray(seq(1,0.1,-0.1))
source('~/Documents/R/color.palette.R')
pal <- color.palette(c("red","yellow","white","cyan","blue"), c(20,5,5,20))
cm=pal(10)
cm[5:6]="white"

annR<-ann<-array(0,dim=c(length(lats),length(lons),length(years)))
for(i in 1:length(lons))
  for(j in 1:length(lats))
    for(k in 1:length(years))
  {
    I=which(ncep1[,4]>=floor(lons[i]) & ncep1[,4]<ceiling(lons[i]) & ncep1[,5]>=floor(lats[j]) & ncep1[,5]<ceiling(lats[j]) & yy==years[k] )
    ann[j,i,k]=length(I)
    annR[j,i,k]=sum(ncep1[I,8]) ## of fixes w/ Morwenna rain on subs day
  }

###Re-do at lower resolution. 
setwd('~/Documents/ECLs')
ncep1<-read.csv('CSV/ECLfixes_unsw_ncep1_250.csv')
ncep1=ncep1[ncep1$CV>=0.25,]
source('~/Documents/R/ECLextract.R')
lats=seq(-40,-24,2)
lons=seq(150,160,2)
cool=ECLcount_seas(ncep1,lats,lons,c(5,10))
warm=ECLcount_seas(ncep1,lats,lons,c(11,4))
all=ECLcount_seas(ncep1,lats,lons,c(1,12))

##Plotty plot plot
cols=gray(seq(1,0.1,-0.1))
colsa=gray(seq(1,0,-0.09))

early=apply(all[,,1:30],c(1,2),mean)
late=apply(all[,,31:62],c(1,2),mean)*30/33
tiff(file=paste("Locs_NCEP1_ann_5079.tiff",sep=""), height=600, width=500)
early[early>6]=6
image.plot(lons,lats,t(early),breaks=seq(0,6,0.5),xlab="",ylab="",col=colsa,zlim=c(0,6))
contour(Useful$x,Useful$y,mask,add=T,drawlabels=F)
dev.off()
tiff(file=paste("Locs_NCEP1_ann_8012.tiff",sep=""), height=600, width=500)
late[late>6]=6
image.plot(lons,lats,t(late),breaks=seq(0,6,0.5),xlab="",ylab="",col=colsa,zlim=c(0,6))
contour(Useful$x,Useful$y,mask,add=T,drawlabels=F)
dev.off()

col2=pal(16)
col2[8:9]="white"
aa=late-early
aa[aa<=(-2)]=-1.95
aa[aa>2]=2
tiff(file=paste("Locs_NCEP1_ann_change.tiff",sep=""), height=600, width=500)
image.plot(lons,lats,t(late-early),breaks=seq(-2,2,0.25),xlab="",ylab="",col=col2,zlim=c(-2,2))
contour(Useful$x,Useful$y,mask,add=T,drawlabels=F)
dev.off()

##ENSO plotty plot plot

read.csv('~/Documents/Timeseries/state.csv',sep=";")->enso

I=which(enso[,2]=="E")
E=apply(cool[,,I],c(1,2),mean)
I=which(enso[,2]=="N")
N=apply(cool[,,I],c(1,2),mean)
I=which(enso[,2]=="L")
L=apply(cool[,,I],c(1,2),mean)

cols=gray(seq(1,0,-0.1))
image.plot(lons,lats,t(E),breaks=seq(0,5.5,0.5),xlab="",ylab="",col=cols,zlim=c(0,5.5))
contour(Useful$x,Useful$y,mask,add=T,drawlabels=F)

I=which(enso[1:62,2]=="E")
E=apply(warm[,,I],c(1,2),mean)
I=which(enso[1:62,2]=="N")
N=apply(warm[,,I],c(1,2),mean)
I=which(enso[1:62,2]=="L")
L=apply(warm[,,I],c(1,2),mean)

cols=gray(seq(1,0.1,-0.1))
image.plot(lons,lats,t(N),breaks=seq(0,3,0.3),xlab="",ylab="",col=cols,zlim=c(0,3))
contour(Useful$x,Useful$y,mask,add=T,drawlabels=F)



