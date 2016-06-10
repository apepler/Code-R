###Figure 3 - Uwind

rm(list=ls())
library(RNetCDF)
setwd("~/Documents/Data/CMIP5")
read.csv("Modellist.csv")->list
list=list[c(-24,-14),]
load("CMIP5corrs.RData")
Ucomp2=cbind(Ucomp[,1:2],apply(Ucomp[,3:39],1,mean,na.rm=T))

tiff(file="~/Documents/Data/CMIP5/Figure3.tiff", height=500, width=800)
mnames=c("Jan","Feb","Mar","Apr","May","Jun","Jul","Aug","Sep","Oct","Nov","Dec")
plot(seq(1,12),Ucomp2[,3],type="l",col="black",xaxt="n", xlab="", ylab="",lwd=2,ylim=c(-5,10),axes=F)
axis(1, at = seq(1:12), labels = mnames,cex.axis=1.5) 
axis(2, at = seq(-5,10,5),cex.axis=1.5)
title(main="Average Uwind (m/s)",cex.main=2)
legend("topleft",legend=c("NCEP","CMIP5 Average"),col=c("blue","red"),lwd=2,bty="n",cex=1.5)
for(i in 1:37) lines(seq(1:12),Ucomp[,i+2],col="gray")
lines(seq(1,12),Ucomp2[,2],col="blue",lwd=2)
lines(seq(1,12),Ucomp2[,3],col="red",lwd=2)
abline(h=0,col="black",lwd=2,lty=2)
dev.off()

#### Figures 5 and 10 - scatter plots

data=read.csv("/home/nfs/z3478332/Documents/Data/CMIP5/25-35/newcomp_v3.csv")
##Only want ONE AWAP version - currently have 2.5 deg (1), 0.5 deg (2), 0.05 deg (3) and 0.05 deg landmask (4)
##Corrs are basically the same, except for 2.5 SW slightly weaker (by ~0.05)
##Let's stick with 2.5 for now
data=data[-c(1,2,3),]
tiff(file="~/Documents/Data/CMIP5/Figure5.tiff",height=600,width=600,pointsize=20)
a=6
b=5
plot(data[2:38,a],data[2:38,b],type="p",col="grey",xlab="East Coast",ylab="South Coast",xlim=c(-0.7,0.2),ylim=c(-0.7,0.7))
abline(h=0,col="grey")
abline(v=0,col="grey")
##Actually, want: Black for AWAP, grey for both right, red for only coast right, unfilled red where neither
I=which(data[,a]<=(-0.265) & data[,b]<0.265)
points(data[I,a],data[I,b],col="darkgrey",pch=16)
I=which(data[,a]<=(-0.265) & data[,b]>=0.265)
points(data[I[2:length(I)],a],data[I[2:length(I)],b],col="darkgrey",pch=18,cex=2)
points(data[1,a],data[1,b],col="black",pch=4,cex=2,lwd=4)
points(median(data[2:38,a]),median(data[2:38,b]),col="black",pch=1,cex=2,lwd=4)
dev.off()

data=read.csv("/home/nfs/z3478332/Documents/Data/CMIP5/25-35/newcomp_v3.csv")
##Only want ONE AWAP version - currently have 2.5 deg (1), 0.5 deg (2), 0.05 deg (3) and 0.05 deg landmask (4)
##Corrs are basically the same, except for 2.5 SW slightly weaker (by ~0.05)
##Let's stick with 2.5 for now
data=data[-c(1,2,3),]
cmip=apply(data[2:38,5:10],2,median)

data1=read.csv("/home/nfs/z3478332/Documents/Data/CMIP5/CMIP3/20c3m/newcomp_wrf50_postregrid.csv")
data1=data1[-c(1,2,3),]
data2=read.csv("/home/nfs/z3478332/Documents/Data/CMIP5/CMIP3/20c3m/newcomp_wrf10_postregrid.csv")
data2=data2[-c(1,2,3),]
a=11
b=10
fname="~/Documents/Data/CMIP5/CMIP3/20c3m/WRF10/Figure10.tiff"
tiff(file=fname,height=600,width=600,pointsize=20)
plot(data2[2:5,a],data2[2:5,b],type="p",col="gray50",pch=19,xlab="East Coast",ylab="South Coast",xlim=c(-0.7,0.2),ylim=c(-0.7,0.7))
abline(h=0,col="grey")
abline(v=0,col="grey")

points(data1[6:17,a],data1[6:17,b],col="gray65",pch=17)
points(data2[6:17,a],data2[6:17,b],col="gray80",pch=15)
points(data1[1,a],data1[1,b],col="black",pch=4,cex=2,lwd=4)
points(median(data1[2:5,a]),median(data1[2:5,b]),col="gray50",pch=1,cex=2,lwd=4)
points(median(data1[6:17,a]),median(data1[6:17,b]),col="gray65",pch=2,cex=2,lwd=4)
points(median(data2[6:17,a]),median(data2[6:17,b]),col="gray80",pch=0,cex=2,lwd=4)
points(cmip[2],cmip[1],col="black",pch=1,cex=2,lwd=4)
legend("bottomright",c("AWAP","CMIP5","CMIP3","WRF50","WRF10"),
       pch=c(4,1,1,2,0),col=c("black","black","gray50","gray65","gray80"),pt.lwd=2)
dev.off()

##### Figure 9

rm(list=ls())
library(fields)
library("R.matlab")
readMat('~/Documents/Data/Useful.mat')->Useful
mask<-t(Useful$mask)
mask[is.na(mask)]=0
source('~/R/x86_64-pc-linux-gnu-library/2.15/gplots/R/rich.colors.R')
source('~/Documents/R/color.palette.R')
pal <- color.palette(c("red","yellow","white","cyan","blue"), c(20,5,5,20))
cm=pal(14)
cm[7:8]="white"
setwd('/home/nfs/z3478332/Documents/Data/CMIP5/CMIP3/')
breaklist=c(seq(-0.7,-0.1,0.1),0,seq(0.1,0.7,0.1))
ColorBar <- function(brks,cols,vert=T,subsampleg=2)
{
  par(mar = c(3, 1, 3, 3), mgp = c(1, 1, 0), las = 1, cex = 1.2)
  image(1, c(1:length(cols)), t(c(1:length(cols))), axes = FALSE, col = cols, 
        xlab = '', ylab = '')
  box()
  axis(4, at = seq(0.5, length(brks) - 0.5, subsampleg), tick = TRUE, 
       labels = brks[seq(1, length(brks), subsampleg)])
}


load('WRFR10corrs.RData')
cool10=apply(coolA,c(1,2),median)
lat10=seq(-45,-20,0.1)
lon10=seq(130,160,0.1)
sigmask10=which(abs(cool10)>=0.45,arr.ind=T)
I=which(lon10[sigmask10[,1]]>155 | lon10[sigmask10[,1]]<135)
if(length(I)>1) sigmask10=sigmask10[-I,]
I=which(lat10[sigmask10[,2]]>(-23) | lat10[sigmask10[,2]]<(-39))
if(length(I)>1) sigmask10=sigmask10[-I,]
I=which(lon10[sigmask10[,1]] %% 0.5 == 0 & lat10[sigmask10[,2]] %% 0.5 == 0)
if(length(I)>1) sigmask10=sigmask10[I,]
cool10[cool10>=0.69]=0.69
cool10[cool10<=(-0.69)]=-0.69

load('WRF50corrs.RData')
cool50=apply(coolA,c(1,2),median)
lat50=seq(-50,0,0.5)
lon50=seq(105,180,0.5)
sigmask50=which(abs(cool50)>=0.45,arr.ind=T)
I=which(lon50[sigmask50[,1]]>155 | lon50[sigmask50[,1]]<135)
if(length(I)>1) sigmask50=sigmask50[-I,]
I=which(lat50[sigmask50[,2]]>(-23) | lat50[sigmask50[,2]]<(-39))
if(length(I)>1) sigmask50=sigmask50[-I,]
cool50[cool50>=0.69]=0.69
cool50[cool50<=(-0.69)]=-0.69

pdf(file="Figure9.pdf", height=7, width=18, pointsize=20)
layout(cbind(1,2,3), width=c(1.1,1,0.3),height=1)
par(mar=c(5,5,1,1)+0.1)
image(lon50,lat50,cool50,xlab="",breaks=seq(-0.7,0.7,0.1),ylab="",col=cm,zlim=c(-0.7,0.7)  ,xlim=c(135,155),ylim=c(-39,-23),cex.axis=2)
contour(Useful$x,Useful$y,mask,add=T,drawlabels=F)
points(lon50[sigmask50[,1]],lat50[sigmask50[,2]],col="black",pch=19,cex=0.75,xlim=c(135,155),ylim=c(-39,-23))
par(mar=c(5,1.5,1,1)+0.1)
image(lon10,lat10,cool10,xlab="",breaks=seq(-0.7,0.7,0.1),ylab="",col=cm,zlim=c(-0.7,0.7)  ,xlim=c(135,155),ylim=c(-39,-23),axes=F)
axis(1,cex.axis=2)
contour(Useful$x,Useful$y,mask,add=T,drawlabels=F)
points(lon10[sigmask10[,1]],lat10[sigmask10[,2]],col="black",pch=19,cex=0.75,xlim=c(135,155),ylim=c(-39,-23))
ColorBar(brks = breaklist, cols = cm)
dev.off()

############ Figure 8 

load('CMIP3corrs.RData')
lat=seq(-90,90,2.5)
lon=seq(0,357.5,2.5)
cool=apply(coolC,c(1,2),median)
sigmask=which(abs(cool)>=0.265,arr.ind=T)
I=which(lon[sigmask[,1]]>160 | lon[sigmask[,1]]<110)
if(length(I)>1) sigmask=sigmask[-I,]
I=which(lat[sigmask[,2]]>(-10) | lat[sigmask[,2]]<(-45))
if(length(I)>1) sigmask=sigmask[-I,]
cool[cool>=0.69]=0.69
cool[cool<=(-0.69)]=-0.69

##Change sigmask for WRF50 so not quite as annoying.
sigmask50=which(abs(cool50)>=0.45,arr.ind=T)
I=which(lon50[sigmask50[,1]]>110 | lon50[sigmask50[,1]]<160)
if(length(I)>1) sigmask=sigmask50[-I,]
I=which(lat50[sigmask50[,2]]>(-10) | lat50[sigmask50[,2]]<(-45))
if(length(I)>1) sigmask50=sigmask50[-I,]
I=which(lon50[sigmask50[,1]] %% 1 == 0 & lat50[sigmask50[,2]] %% 1 == 0)
if(length(I)>1) sigmask50=sigmask50[I,]

pdf(file="Figure8.pdf", height=7, width=18, pointsize=20)
layout(cbind(1,2,3), width=c(1.1,1,0.3),height=1)
par(mar=c(5,5,1,1)+0.1)
image(lon,lat,cool,xlab="",breaks=seq(-0.7,0.7,0.1),ylab="",col=cm,zlim=c(-0.7,0.7),xlim=c(110,160),ylim=c(-45,-10),cex.axis=2)
contour(Useful$x,Useful$y,mask,add=T,drawlabels=F)
points(lon[sigmask[,1]],lat[sigmask[,2]],col="black",pch=19,cex=0.75,xlim=c(110,160),ylim=c(-45,-10))
par(mar=c(5,1.5,1,1)+0.1)
image(lon50,lat50,cool50,xlab="",breaks=seq(-0.7,0.7,0.1),ylab="",col=cm,zlim=c(-0.7,0.7),xlim=c(110,160),ylim=c(-45,-10),axes=F)
axis(1,cex.axis=2)
contour(Useful$x,Useful$y,mask,add=T,drawlabels=F)
points(lon50[sigmask50[,1]],lat50[sigmask50[,2]],col="black",pch=19,cex=0.5,xlim=c(110,160),ylim=c(-45,-10))
ColorBar(brks = breaklist, cols = cm)
dev.off()

####Figure 7

load('NCEP50corrs.RData')
cool50=apply(coolA,c(1,2),median)
cool50[cool50>=0.69]=0.69
cool50[cool50<=(-0.69)]=-0.69
sigmask50=which(abs(cool50)>=0.265,arr.ind=T)
I=which(lon50[sigmask50[,1]]>110 | lon50[sigmask50[,1]]<160)
if(length(I)>1) sigmask=sigmask50[-I,]
I=which(lat50[sigmask50[,2]]>(-10) | lat50[sigmask50[,2]]<(-45))
if(length(I)>1) sigmask50=sigmask50[-I,]
I=which(lon50[sigmask50[,1]] %% 1 == 0 & lat50[sigmask50[,2]] %% 1 == 0)
if(length(I)>1) sigmask50=sigmask50[I,]

pdf(file="Figure7.pdf", height=7, width=11, pointsize=20)
layout(cbind(1,2), width=c(1.1,0.3),height=1)
par(mar=c(3,3,1,1)+0.1)
image(lon50,lat50,cool50,xlab="",breaks=seq(-0.7,0.7,0.1),ylab="",col=cm,zlim=c(-0.7,0.7),xlim=c(110,160),ylim=c(-45,-10),cex.axis=1.2)
contour(Useful$x,Useful$y,mask,add=T,drawlabels=F)
points(lon50[sigmask50[,1]],lat50[sigmask50[,2]],col="black",pch=19,cex=0.3,xlim=c(110,160),ylim=c(-45,-10))
ColorBar(brks = breaklist, cols = cm)
dev.off()

##### Figure 6

setwd('/home/nfs/z3478332/Documents/Data/CMIP5/')
load("CMIP5corrs.RData")
lat=seq(-90,90,2.5)
lon=seq(0,357.5,2.5)
good=c(22,24,26,27,28,34,35)
cool=apply(coolC[,,good],c(1,2),median)
sigmask=which(abs(cool)>=0.265,arr.ind=T)
I=which(lon[sigmask[,1]]>160 | lon[sigmask[,1]]<110)
if(length(I)>1) sigmask=sigmask[-I,]
I=which(lat[sigmask[,2]]>(-10) | lat[sigmask[,2]]<(-45))
if(length(I)>1) sigmask=sigmask[-I,]
cool[cool>=0.69]=0.69
cool[cool<=(-0.69)]=-0.69

pdf(file="Figure6.pdf", height=7, width=11, pointsize=20)
layout(cbind(1,2), width=c(1.1,0.3),height=1)
par(mar=c(3,3,1,1)+0.1)
image(lon,lat,cool,xlab="",breaks=seq(-0.7,0.7,0.1),ylab="",col=cm,zlim=c(-0.7,0.7),xlim=c(110,160),ylim=c(-45,-10),cex.axis=1.2)
contour(Useful$x,Useful$y,mask,add=T,drawlabels=F)
points(lon[sigmask[,1]],lat[sigmask[,2]],col="black",pch=19,cex=0.5,xlim=c(110,160),ylim=c(-45,-10))
ColorBar(brks = breaklist, cols = cm)
dev.off()

#### Figure 4
library(RNetCDF)
setwd('/home/nfs/z3478332/Documents/Data/CMIP5/')
load("CMIP5corrs.RData")
lat=seq(-90,90,2.5)
lon=seq(0,357.5,2.5)

coolC=apply(coolC,c(1,2),median)
sigmaskC=which(abs(coolC)>=0.265,arr.ind=T)
I=which(lon[sigmaskC[,1]]>160 | lon[sigmaskC[,1]]<110)
if(length(I)>1) sigmaskC=sigmaskC[-I,]
I=which(lat[sigmaskC[,2]]>(-10) | lat[sigmaskC[,2]]<(-45))
if(length(I)>1) sigmaskC=sigmaskC[-I,]
coolC[coolC>=0.69]=0.69
coolC[coolC<=(-0.69)]=-0.69

warmC=apply(warmC,c(1,2),median)
sigmaskW=which(abs(warmC)>=0.265,arr.ind=T)
I=which(lon[sigmaskW[,1]]>160 | lon[sigmaskW[,1]]<110)
if(length(I)>1) sigmaskW=sigmaskW[-I,]
I=which(lat[sigmaskW[,2]]>(-10) | lat[sigmaskW[,2]]<(-45))
if(length(I)>1) sigmaskW=sigmaskW[-I,]
warmC[warmC>=0.69]=0.69
warmC[warmC<=(-0.69)]=-0.69

load('~/Documents/Data/AWAPrain.RData')
yearsA=seq(1900,2012)
I=which(yearsA>=1950 & yearsA<=2005)
cool=cool[,,I]
yearsB=seq(1901,2012)
I=which(yearsB>=1950 & yearsB<=2005) #So has the year ending in 1950 - year ending in 2005
warm=warm[,,I] 

f1=open.nc('/srv/ccrc/data34/z3478332/NCEP/Monthly/uwnd.mon.mean.nc')
timeN=var.get.nc(f1,'time')
timeN=as.Date((timeN/24-2),origin="0001-01-01")
timeN=cbind(as.numeric(format(timeN,"%Y")),as.numeric(format(timeN,"%m")))
##If using northern region
NCEP=apply(var.get.nc(f1,'uwnd',c(61,47,3,1),c(1,5,1,length(timeN[,1])),unpack=T),2,mean,na.rm=T)
years=cbind(seq(1950,2005),matrix(0,56,2))
for(i in 1:length(years[,1]))
{
  I=which(timeN[,1]==years[i,1] & timeN[,2]>=5 & timeN[,2]<=10)
  years[i,2]=mean(NCEP[I])
  I=which((timeN[,1]==years[i,1] & timeN[,2]<=4) | (timeN[,1]==years[i,1]-1 & timeN[,2]>=11))
  years[i,3]=mean(NCEP[I])
}

coolA<-warmA<-matrix(0,691,886)
for(i in 1:691)
  for(j in 1:886)
  {
    coolA[i,j]=cor(cool[i,j,],years[,2],use="pairwise.complete.obs")
    warmA[i,j]=cor(warm[i,j,],years[,3],use="pairwise.complete.obs")
  }
coolA[coolA>=0.69]=0.69
coolA[coolA<=(-0.69)]=-0.69
warmA[warmA>=0.69]=0.69
warmA[warmA<=(-0.69)]=-0.69

reglats=c(rep(-32.5,4),rep(-35,4),rep(-37.5,3))
reglons=c(seq(135,142.5,2.5),seq(137.5,145,2.5),seq(140,145,2.5))
reglats2=seq(-25,-32.5,-2.5)

#coolA[abs(coolA)<=0.265]=NaN
#warmA[abs(warmA)<=0.265]=NaN

pdf(file="Figure4_all2.pdf", height=14, width=18,pointsize=20)
layout(cbind(c(1,3),c(2,4),c(5,6)), width=c(1.1,1,0.3),height=c(1,1))

par(mar=c(1.5,5,5,1)+0.1)
image(as.vector(Useful$x),as.vector(Useful$y),t(warmA*Useful$mask),xlab="",
      breaks=seq(-0.7,0.7,0.1),ylab="",col=cm,zlim=c(-0.7,0.7),xlim=c(110,160),ylim=c(-45,-10),axes=F)
axis(2,cex.axis=2)
text(115,-15,"a)",cex=2)
contour(Useful$x,Useful$y,mask,add=T,drawlabels=F)

par(mar=c(1.5,1.5,5,1)+0.1)
image(as.vector(Useful$x),as.vector(Useful$y),t(coolA*Useful$mask),xlab="",
      breaks=seq(-0.7,0.7,0.1),ylab="",col=cm,zlim=c(-0.7,0.7),xlim=c(110,160),ylim=c(-45,-10),axes=F)
text(115,-15,"b)",cex=2)
contour(Useful$x,Useful$y,mask,add=T,drawlabels=F)
points(reglons,reglats,pch=4,col="gray",lwd=4,cex=2)
polygon(x=c(131.25,141.25,146.75,146.75,138.75,131.25),y=c(-31.25,-31.25,-36.25,-38.75,-38.75,-31.25),lwd=4,border="gray")
points(rep(152.5,4),reglats2,pch=4,col="black",lwd=4,cex=2)
rect(151.25,-33.75,153.75,-23.75,lwd=4)

par(mar=c(5,5,1,1)+0.1)
image(lon,lat,warmC,xlab="",breaks=seq(-0.7,0.7,0.1),ylab="",col=cm,zlim=c(-0.7,0.7)  ,xlim=c(110,160),ylim=c(-45,-10),cex.axis=2)
text(115,-15,"c)",cex=2)
contour(Useful$x,Useful$y,mask,add=T,drawlabels=F)
points(lon[sigmaskW[,1]],lat[sigmaskW[,2]],col="black",pch=19,cex=0.75,xlim=c(110,160),ylim=c(-45,-10))

par(mar=c(5,1.5,1,1)+0.1)
image(lon,lat,coolC,xlab="",breaks=seq(-0.7,0.7,0.1),ylab="",col=cm,zlim=c(-0.7,0.7),xlim=c(110,160),ylim=c(-45,-10),axes=F)
text(115,-15,"d)",cex=2)
axis(1,cex.axis=2)
contour(Useful$x,Useful$y,mask,add=T,drawlabels=F)
points(lon[sigmaskC[,1]],lat[sigmaskC[,2]],col="black",pch=19,cex=0.75,xlim=c(110,160),ylim=c(-45,-10))
ColorBar(brks = breaklist, cols = cm)
ColorBar(brks = breaklist, cols = cm)
dev.off()

###Figure 2 - is a problem.

rm(list=ls())
library(RNetCDF)
setwd("~/Documents/Data/CMIP5")
read.csv("Modellist.csv")->list
list=list[c(-24,-14),]
source('~/Documents/R/CMIP5_play.R')

library(fields)
source('~/R/x86_64-pc-linux-gnu-library/2.15/gplots/R/rich.colors.R')
source('~/Documents/R/color.palette.R')
pal <- color.palette(c("red","yellow","white","cyan","blue"), c(20,5,5,20))
cm=pal(14)
ColorBar <- function(brks,cols,vert=T,subsampleg=1)
{
  par(mar = c(1, 1, 1, 3), mgp = c(1, 1, 0), las = 1, cex = 1.2)
  image(1, c(1:length(cols)), t(c(1:length(cols))), axes = FALSE, col = cols, 
        xlab = '', ylab = '')
  box()
  axis(4, at = seq(1.5, length(brks) - 1.5, subsampleg), tick = TRUE, 
       labels = brks[seq(2, length(brks)-1, subsampleg)])
}

cool<-warm<-array(0,c(144,73,37))
for(i in 1:37)
{
  a=CMIP_rain_v2(as.character(list[i,1]),as.character(list[i,2]),as.character(list[i,3]),"N")
  cool[,,i]=unlist(a[2])
  warm[,,i]=unlist(a[3])
}
coolC=apply(cool,c(1,2),mean,na.rm=T)
warmC=apply(warm,c(1,2),mean,na.rm=T)

load('~/Documents/Data/AWAPrain.RData')
yearsA=seq(1900,2012)
I=which(yearsA>=1950 & yearsA<=2005)
coolA=apply(cool[,,I],c(1,2),mean)
yearsB=seq(1901,2012)
I=which(yearsB>=1950 & yearsB<=2005) #So has the year ending in 1950 - year ending in 2005
warmA=apply(warm[,,I],c(1,2),mean) #So has the year ending in 1950 - year ending in 2005
library(fields)
a=interp.surface.grid(list(x=Useful$x,y=Useful$y,z=t(coolA)),list(x=lon,y=lat))
coolA1=a$z
a=interp.surface.grid(list(x=Useful$x,y=Useful$y,z=t(warmA)),list(x=lon,y=lat))
warmA1=a$z

##Now need to re-do my stuff on creating a semi-mask for WRF

library("R.matlab")
readMat('~/Documents/Data/Useful.mat')->Useful
mask<-t(Useful$mask)
mask[is.na(mask)]=0

maskC<-maskCa<-matrix(0,144,73)
for(i in 2:143)
  for(j in 2:72)
  {
    I=which(Useful$x>lon[i]-1.25 & Useful$x<lon[i]+1.25)
    J=which(Useful$y>lat[j]-1.25 & Useful$y<lat[j]+1.25)
    
    if(length(I)>0 & length(J)>0) maskC[i,j]=mean(mask[I,J])
  }
maskCa[maskC>=0.25]=1
maskCb=maskCa
maskCb[maskCa<1]=NaN

diffC=coolC-coolA1
diffW=warmC-warmA1
pcC=(diffC/coolA1)*100
pcW=(diffW/warmA1)*100

## Just anomaly together

breaklist=c(-1000,-200,-150,-100,-50,-25,-10,0,10,25,50,100,150,200,1000)
pdf(file="Figure2_anoms.pdf", height=7, width=18, pointsize=20)
layout(cbind(1,2,3), width=c(1.1,1,0.3),height=1)
par(mar=c(5,5,1,1)+0.1)
image(lon,lat,diffW*maskCb,xlab="",breaks=breaklist,ylab="",col=cm,zlim=c(-1000,1000),xlim=c(110,160),ylim=c(-45,-10),axes=F)
axis(1,cex.axis=2)
axis(2,cex.axis=2)
contour(Useful$x,Useful$y,mask,add=T,drawlabels=F)
par(mar=c(5,1.5,1,1)+0.1)
image(lon,lat,diffC*maskCb,xlab="",breaks=breaklist,ylab="",col=cm,zlim=c(-1000,1000),xlim=c(110,160),ylim=c(-45,-10),axes=F)
axis(1,cex.axis=2)
contour(Useful$x,Useful$y,mask,add=T,drawlabels=F)
ColorBar(brks = breaklist, cols = cm)
dev.off()


## Just PC change together
ColorBar2 <- function(brks,cols,vert=T,subsampleg=1,brklabs=brks[seq(2, length(brks)-1, subsampleg)])
{
  par(mar = c(1, 1, 1, 4), mgp = c(1, 1, 0), las = 1, cex = 1.2)
  image(1, c(1:length(cols)), t(c(1:length(cols))), axes = FALSE, col = cols, 
        xlab = '', ylab = '')
  box()
  axis(4, at = seq(1.5, length(brks) - 1.5, subsampleg), tick = TRUE, 
       labels = brklabs)
}



breaklist=c(-1000,-100,-80,-60,-40,-20,-10,0,10,20,40,60,80,100,1000)
pdf(file="Figure2_PCchange_v2.pdf", height=7, width=18, pointsize=20)
layout(cbind(1,2,3), width=c(1.1,1,0.3),height=1)
par(mar=c(5,5,1,1)+0.1)
image(lon,lat,pcW*maskCb,xlab="",breaks=breaklist,ylab="",col=cm,zlim=c(-1000,1000),xlim=c(110,160),ylim=c(-45,-10),axes=F)
axis(1,cex.axis=2)
axis(2,cex.axis=2)
contour(Useful$x,Useful$y,mask,add=T,drawlabels=F)
par(mar=c(5,1.5,1,1)+0.1)
image(lon,lat,pcC*maskCb,xlab="",breaks=breaklist,ylab="",col=cm,zlim=c(-1000,1000),xlim=c(110,160),ylim=c(-45,-10),axes=F)
axis(1,cex.axis=2)
contour(Useful$x,Useful$y,mask,add=T,drawlabels=F)
ColorBar2(brks = breaklist, cols = cm,brklabs=c("-100%","-80%","-60%","-40%","-20%","-10%","","+10%","+20%","+40%","+60%","+80%","+100%"))
dev.off()


## Weird combo

sigmaskc1=which(abs(pcC*maskCb)>10 & abs(pcC*maskCb)<=25 ,arr.ind=T)
sigmaskc2=which(abs(pcC*maskCb)>25 & abs(pcC*maskCb)<=100 ,arr.ind=T)
sigmaskc3=which(abs(pcC*maskCb)>100 ,arr.ind=T)
sigmaskw1=which(abs(pcW*maskCb)>10 & abs(pcW*maskCb)<=25 ,arr.ind=T)
sigmaskw2=which(abs(pcW*maskCb)>25 & abs(pcW*maskCb)<=100 ,arr.ind=T)
sigmaskw3=which(abs(pcW*maskCb)>100 ,arr.ind=T)

breaklist=c(-1000,-200,-150,-100,-50,-25,-10,0,10,25,50,100,150,200,1000)
pdf(file="Figure2_combo.pdf", height=7, width=18, pointsize=20)
layout(cbind(1,2,3), width=c(1.1,1,0.3),height=1)
par(mar=c(5,5,1,1)+0.1)
image(lon,lat,diffW*maskCb,xlab="",breaks=breaklist,ylab="",col=cm,zlim=c(-1000,1000),xlim=c(110,160),ylim=c(-45,-10),cex.axis=2)
axis(1,cex.axis=2)
axis(2,cex.axis=2)
points(lon[sigmaskw1[,1]],lat[sigmaskw1[,2]],col="white",pch=19,cex=0.5,xlim=c(110,160),ylim=c(-45,-10))
points(lon[sigmaskw2[,1]],lat[sigmaskw2[,2]],col="grey",pch=19,cex=0.5,xlim=c(110,160),ylim=c(-45,-10))
points(lon[sigmaskw3[,1]],lat[sigmaskw3[,2]],col="black",pch=19,cex=0.5,xlim=c(110,160),ylim=c(-45,-10))
contour(Useful$x,Useful$y,mask,add=T,drawlabels=F)

par(mar=c(5,1.5,1,1)+0.1)
image(lon,lat,diffC*maskCb,xlab="",breaks=breaklist,ylab="",col=cm,zlim=c(-1000,1000),xlim=c(110,160),ylim=c(-45,-10),axes=F)
axis(1,cex.axis=2)
points(lon[sigmaskc1[,1]],lat[sigmaskc1[,2]],col="white",pch=19,cex=0.5,xlim=c(110,160),ylim=c(-45,-10))
points(lon[sigmaskc2[,1]],lat[sigmaskc2[,2]],col="grey",pch=19,cex=0.5,xlim=c(110,160),ylim=c(-45,-10))
points(lon[sigmaskc3[,1]],lat[sigmaskc3[,2]],col="black",pch=19,cex=0.5,xlim=c(110,160),ylim=c(-45,-10))
contour(Useful$x,Useful$y,mask,add=T,drawlabels=F)
ColorBar(brks = breaklist, cols = cm)
dev.off()

###Normalised?

load('~/Documents/Data/AWAPrain.RData')
yearsA=seq(1900,2012)
I=which(yearsA>=1950 & yearsA<=2005)
yearsB=seq(1901,2012)
J=which(yearsB>=1950 & yearsB<=2005) #So has the year ending in 1950 - year ending in 2005

cool2<-warm2<-array(NaN,c(144,73,length(I)))

for(i in 1:length(I))
{
  a=interp.surface.grid(list(x=Useful$x,y=Useful$y,z=t(cool[,,I[i]])),list(x=lon,y=lat))
  cool2[,,i]=a$z
  a=interp.surface.grid(list(x=Useful$x,y=Useful$y,z=t(warm[,,J[i]])),list(x=lon,y=lat))
  warm2[,,i]=a$z
}

coolA2=apply(cool2,c(1,2),mean)
coolAsd=apply(cool2,c(1,2),sd)
warmA2=apply(warm2,c(1,2),mean)
warmAsd=apply(warm2,c(1,2),sd)

diffCnorm=(coolC-coolA2)/coolAsd
diffWnorm=(warmC-warmA2)/warmAsd

breaklist=c(-5,-3,-2,-1.5,-1,-0.5,0,0.5,1,1.5,2,3,5)
pdf(file="Figure2_normchange_v2.pdf", height=7, width=18, pointsize=20)
layout(cbind(1,2,3), width=c(1.1,1,0.3),height=1)
par(mar=c(5,5,1,1)+0.1)
image(lon,lat,diffWnorm*maskCb,xlab="",breaks=breaklist,ylab="",col=cm[2:13],zlim=c(-5,5),xlim=c(110,160),ylim=c(-45,-10),axes=F)
axis(1,cex.axis=2)
axis(2,cex.axis=2)
contour(Useful$x,Useful$y,mask,add=T,drawlabels=F)
par(mar=c(5,1.5,1,1)+0.1)
image(lon,lat,diffCnorm*maskCb,xlab="",breaks=breaklist,ylab="",col=cm[2:13],zlim=c(-5,5),xlim=c(110,160),ylim=c(-45,-10),axes=F)
axis(1,cex.axis=2)
contour(Useful$x,Useful$y,mask,add=T,drawlabels=F)
ColorBar(brks = breaklist, cols = cm[2:13])
dev.off()


