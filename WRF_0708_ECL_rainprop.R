rm(list=ls())
setwd('~/Documents/ECLs/WRFruns/0708/')

library(RNetCDF)
library(fields)
library(akima)
library("R.matlab")
readMat('~/Documents/Data/Useful_ECL.mat')->UsefulECL
readMat('~/Documents/Data/Useful.mat')->Useful
mask<-t(Useful$mask)
mask[is.na(mask)]=0
mask2<-t(UsefulECL$mask)
mask2[is.na(mask2)]=0
source('~/R/x86_64-pc-linux-gnu-library/2.15/gplots/R/rich.colors.R')
wrfv=c("R1","R2","R3")

cat=c("rad5_p100","rad5_p240","rad2_p100","rad2_p240")
c=3

prop<-propNT<-array(0,c(144,215,3))
prop2<-propNT2<-array(0,c(886,691,3))

f1=open.nc("/srv/ccrc/data13/z3393020/Analyses/share/geo_em_files/geo_em.d01.narclim.nc")
lat1=var.get.nc(f1,"XLAT_M")
lon1=var.get.nc(f1,"XLONG_M")
latt=as.vector(lat1)
lont=as.vector(lon1)
library(akima)

for(r in 1:3)
{
  a=open.nc(paste("/srv/ccrc/data36/z3478332/WRF/output/ERAI_R",r,"_nudging_default_2007/out/ECLrain_bymonth_0708_",cat[c],".nc",sep=""))
  ECLrain=var.get.nc(a,"ECLrain_loc")
  rain=var.get.nc(a,"allrain")
  prop[,,r]=apply(ECLrain,c(1,2),sum,na.rm=T)/apply(rain,c(1,2),sum,na.rm=T)
  rr=as.vector(prop[,,r])
  rr[which(is.na(rr))]=0
  rr2=interp(lont,latt,rr,Useful$x,Useful$y)
  prop2[,,r]=rr2$z
  
  a=open.nc(paste("/srv/ccrc/data36/z3478332/WRF/output/ERAI_R",r,"_nudging_default_2007_notopo/out/ECLrain_bymonth_0708_",cat[c],".nc",sep=""))
  ECLrain=var.get.nc(a,"ECLrain_loc")
  rain=var.get.nc(a,"allrain")
  propNT[,,r]=apply(ECLrain,c(1,2),sum,na.rm=T)/apply(rain,c(1,2),sum,na.rm=T)
  rr=as.vector(propNT[,,r])
  rr[which(is.na(rr))]=0
  rr2=interp(lont,latt,rr,Useful$x,Useful$y)
  propNT2[,,r]=rr2$z
  }

source('~/Documents/R/filled.contour3.R')
source('~/Documents/R/filled.legend.R')
prop3=apply(prop2,c(1,2),mean)
prop3[prop3==0]=NaN
propNT3=apply(propNT2,c(1,2),mean)
propNT3[propNT3==0]=NaN
ColorBar <- function(brks,cols,vert=T,subsampleg=1,blab=brks)
{
  par(mar = c(3, 1, 3, 3), mgp = c(1, 1, 0), las = 1, cex = 1)
  image(1, c(1:length(cols)), t(c(1:length(cols))), axes = FALSE, col = cols, 
        xlab = '', ylab = '')
  box()
  axis(4, at = seq(1.5, length(brks) - 1.5, subsampleg), tick = TRUE, 
       labels = blab[seq(2, length(brks)-1, subsampleg)])
}

pdf(file="Rainprop_6hrECL_mean.pdf",width=6,height=4)
layout(cbind(1,2,3),width=c(1,1,0.4))
par(mar=c(2,2,4,2))
filled.contour3(Useful$x,Useful$y,prop3,levels=seq(0,0.4,0.05),col=tim.colors(8),xlim=c(145,156),ylim=c(-40,-20),main="Control")
contour(Useful$x,Useful$y,mask,add=T,drawlabels=F)
contour(Useful$x,Useful$y,mask2,add=T,drawlabels=F)
filled.contour3(Useful$x,Useful$y,propNT3,levels=seq(0,0.4,0.05),col=tim.colors(8),xlim=c(145,156),ylim=c(-40,-20),main="NoTopo")
contour(Useful$x,Useful$y,mask,add=T,drawlabels=F)
contour(Useful$x,Useful$y,mask2,add=T,drawlabels=F)
ColorBar(brks=seq(0,0.4,0.05),col=tim.colors(8),blab=paste(seq(0,40,5),"%",sep=""))
dev.off()

pdf(file="Rainprop_6hrECL_control.pdf",width=9,height=4)
layout(cbind(1,2,3,4),width=c(1,1,1,0.4))
par(mar=c(2,2,4,2))
for(r in 1:3)
{
filled.contour3(Useful$x,Useful$y,prop2[,,r],levels=seq(0,0.4,0.05),col=tim.colors(8),xlim=c(145,156),ylim=c(-40,-20),main=paste("R",r,sep=""))
contour(Useful$x,Useful$y,mask,add=T,drawlabels=F)
contour(Useful$x,Useful$y,mask2,add=T,drawlabels=F)
}
ColorBar(brks=seq(0,0.4,0.05),col=tim.colors(8),blab=paste(seq(0,40,5),"%",sep=""))
dev.off()

########## Repeat for warmcold

propW<-propC<-array(0,c(144,215,6))
propW2<-propC2<-array(0,c(886,691,6))

for(r in 1:3)
{
  a=open.nc(paste("/srv/ccrc/data36/z3478332/WRF/output/ERAI_R",r,"_nudging_default_2007/out/ECLrain_bymonth_0708_",cat[c],".nc",sep=""))
  ECLrain=var.get.nc(a,"ECLrain_loc")
  rain=var.get.nc(a,"allrain")
  propW[,,r]=apply(ECLrain[,,c(1:4,11:16,23:24)],c(1,2),sum,na.rm=T)/apply(rain[,,c(1:4,11:16,23:24)],c(1,2),sum,na.rm=T)
  propC[,,r]=apply(ECLrain[,,c(5:10,17:22)],c(1,2),sum,na.rm=T)/apply(rain[,,c(5:10,17:22)],c(1,2),sum,na.rm=T)

  a=open.nc(paste("/srv/ccrc/data36/z3478332/WRF/output/ERAI_R",r,"_nudging_default_2007_notopo/out/ECLrain_bymonth_0708_",cat[c],".nc",sep=""))
  ECLrain=var.get.nc(a,"ECLrain_loc")
  rain=var.get.nc(a,"allrain")
  propW[,,r+3]=apply(ECLrain[,,c(1:4,11:16,23:24)],c(1,2),sum,na.rm=T)/apply(rain[,,c(1:4,11:16,23:24)],c(1,2),sum,na.rm=T)
  propC[,,r+3]=apply(ECLrain[,,c(5:10,17:22)],c(1,2),sum,na.rm=T)/apply(rain[,,c(5:10,17:22)],c(1,2),sum,na.rm=T)
}

for(r in 1:6)
{
  rr=as.vector(propW[,,r])
  rr[which(is.na(rr))]=0
  rr2=interp(lont,latt,rr,Useful$x,Useful$y)
  propW2[,,r]=rr2$z
  
  rr=as.vector(propC[,,r])
  rr[which(is.na(rr))]=0
  rr2=interp(lont,latt,rr,Useful$x,Useful$y)
  propC2[,,r]=rr2$z
}

prop3=apply(propW2[,,1:3],c(1,2),mean)
prop3[prop3==0]=NaN
propNT3=apply(propW2[,,4:6],c(1,2),mean)
propNT3[propNT3==0]=NaN

pdf(file="Rainprop_6hrECL_mean_warm.pdf",width=6,height=4)
layout(cbind(1,2,3),width=c(1,1,0.4))
par(mar=c(2,2,4,2))
filled.contour3(Useful$x,Useful$y,prop3,levels=seq(0,0.4,0.05),col=tim.colors(8),xlim=c(145,156),ylim=c(-40,-20),main="Control")
contour(Useful$x,Useful$y,mask,add=T,drawlabels=F)
contour(Useful$x,Useful$y,mask2,add=T,drawlabels=F)
filled.contour3(Useful$x,Useful$y,propNT3,levels=seq(0,0.4,0.05),col=tim.colors(8),xlim=c(145,156),ylim=c(-40,-20),main="NoTopo")
contour(Useful$x,Useful$y,mask,add=T,drawlabels=F)
contour(Useful$x,Useful$y,mask2,add=T,drawlabels=F)
ColorBar(brks=seq(0,0.4,0.05),col=tim.colors(8),blab=paste(seq(0,40,5),"%",sep=""))
dev.off()


#### Repeat for SSTs

rm(list=ls())
setwd('~/Documents/ECLs/WRFruns/0708/')

library(RNetCDF)
library(fields)
library(akima)
library("R.matlab")
readMat('~/Documents/GDI/Useful_ECL.mat')->UsefulECL
readMat('~/Documents/GDI/Useful.mat')->Useful
mask<-t(Useful$mask)
mask[is.na(mask)]=0
mask2<-t(UsefulECL$mask)
mask2[is.na(mask2)]=0
source('~/R/x86_64-pc-linux-gnu-library/2.15/gplots/R/rich.colors.R')
wrfv=c("R1","R2","R3")

cat=c("rad5_p100","rad5_p240","rad2_p100","rad2_p240")
c=3

prop<-array(0,c(144,215,3,3))
prop2<-array(0,c(886,691,3,3))

f1=open.nc("/srv/ccrc/data13/z3393020/Analyses/share/geo_em_files/geo_em.d01.narclim.nc")
lat1=var.get.nc(f1,"XLAT_M")
lon1=var.get.nc(f1,"XLONG_M")
latt=as.vector(lat1)
lont=as.vector(lon1)
library(akima)

tlist=c("","_BRAN","_BRAN_noeac","_BRAN_2eac")
dir1=c(36,37,37,36)

for(r in 1:3)
  for(type in 1:3)
{
  a=open.nc(paste("/srv/ccrc/data",dir1[type],"/z3478332/WRF/output/ERAI_R",r,"_nudging_default_2007",tlist[type],"/out/ECLrain_bymonth_0708_",cat[c],".nc",sep=""))
  ECLrain=var.get.nc(a,"ECLrain_loc")
  rain=var.get.nc(a,"allrain")
  prop[,,r,type]=apply(ECLrain,c(1,2),sum,na.rm=T)/apply(rain,c(1,2),sum,na.rm=T)
  rr=as.vector(prop[,,r,type])
  rr[which(is.na(rr))]=0
  rr2=interp(lont,latt,rr,Useful$x,Useful$y)
  prop2[,,r,type]=rr2$z
}

source('~/Documents/R/filled.contour3.R')
source('~/Documents/R/filled.legend.R')
prop3=apply(prop2,c(1,2,4),mean)
prop3[prop3==0]=NaN

ColorBar <- function(brks,cols,vert=T,subsampleg=1,blab=brks)
{
  par(mar = c(3, 1, 3, 3), mgp = c(1, 1, 0), las = 1, cex = 1)
  image(1, c(1:length(cols)), t(c(1:length(cols))), axes = FALSE, col = cols, 
        xlab = '', ylab = '')
  box()
  axis(4, at = seq(1.5, length(brks) - 1.5, subsampleg), tick = TRUE, 
       labels = blab[seq(2, length(brks)-1, subsampleg)])
}

pdf(file="Rainprop_6hrECL_mean_SST.pdf",width=9,height=4)
mlist=c("Control","BRAN","NoEAC")
layout(cbind(1,2,3,4),width=c(1,1,1,0.4))
par(mar=c(2,2,4,2))
for(i in 1:3)
{
filled.contour3(Useful$x,Useful$y,prop3[,,i],levels=seq(0,0.4,0.05),col=tim.colors(8),xlim=c(145,156),ylim=c(-40,-20),main=mlist[i])
contour(Useful$x,Useful$y,mask,add=T,drawlabels=F)
contour(Useful$x,Useful$y,mask2,add=T,drawlabels=F)
}
ColorBar(brks=seq(0,0.4,0.05),col=tim.colors(8),blab=paste(seq(0,40,5),"%",sep=""))
dev.off()

####### What about dividing the change in rainfall into the ECL component & non-ECL component?




