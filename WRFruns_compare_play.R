setwd("~/Documents/ECLs")
library(RNetCDF)
library("R.matlab")
readMat('~/Documents/Data/Useful.mat')->Useful
mask<-t(Useful$mask)
mask[is.na(mask)]=0

trackdir="/home/nfs/z3478332/output/outputUM_wrf_2007-06/"
wrfdir="/srv/ccrc/data36/z3478332/WRF/output/ERAI_R2_nudging_default/out/"
a=open.nc(paste(wrfdir,"WRF_d01_slp_regrid.nc",sep=""))
time=var.get.nc(a,"Times")
lat=var.get.nc(a,"lat0")
lon=var.get.nc(a,"lon0")
slp=var.get.nc(a,"slp0")

list=list(read.csv(paste(trackdir,"ECLfixes_default_fix.csv",sep="")),
          read.csv(paste(trackdir,"ECLfixes_default_50.csv",sep="")),
          read.csv(paste(trackdir,"ECLfixes_default_1deg_fix.csv",sep="")),
          read.csv(paste(trackdir,"ECLfixes_notopo.csv",sep="")))
names=c("Default","Topo/2","1deg Topo","NoTopo")
cm=c("red","blue","darkgreen","purple")
keydates=c(20070609,20070619,20070627)

par(mar=c(4,4,2,2))
for(i in 1:4)
{
  b=list[[i]]
  d2=b$Date+(as.numeric(b$Time)-1)/4
  b$CV[b$Location==0]=NaN
  a=unique(b$ID)
  for(j in 1:length(a)) 
  {
    if(i==1 & j==1)
      plot(d2[b$ID==a[j]],b$CV[b$ID==a[j]],lwd=2,col=cm[i],type="l",
           xlim=c(20070601,20070630),ylim=c(0,4),
           xlab="Date",ylab="Curvature")
    else lines(d2[b$ID==a[j]],b$CV[b$ID==a[j]],lwd=2,col=cm[i])
  }
}
abline(v=keydates,col="gray")
legend("topleft",legend=names,lwd=2,col=cm)

par(mar=c(4,4,2,2))
for(i in 1:4)
{
  b=list[[i]]
  d2=b$Date+(as.numeric(b$Time)-1)/4
  b$MSLP[b$Location==0]=NaN
  a=unique(b$ID)
  for(j in 1:length(a)) 
  {
    if(i==1 & j==1)
      plot(d2[b$ID==a[j]],b$MSLP[b$ID==a[j]],lwd=2,col=cm[i],type="l",
           xlim=c(20070601,20070630),ylim=c(970,1030),
           xlab="Date",ylab="MSLP")
    else lines(d2[b$ID==a[j]],b$MSLP[b$ID==a[j]],lwd=2,col=cm[i])
  }
}
abline(v=keydates,col="gray")
legend("bottomleft",legend=names,lwd=2,col=cm)

pdf(file="WRFruns_Jun07_loc_v2.pdf",width=8,height=3)
layout(cbind(1,2,3))
for(t in 1:length(keydates))
{
  par(mar=c(2,2,2,2))
I=which(time==keydates[t]*100)
contour(Useful$x,Useful$y,mask,drawlabels=F,xlim=c(145,180),ylim=c(-45,-10),main=paste(keydates[t],"00Z"))
contour(lon,lat,slp[,,I],levels=seq(992,1050,4),add=T) #19980808 = 29 (1 + 7*4)
for(i in 1:4)
{
  b=list[[i]]
  I=which(b$Date==keydates[t] & b$Time=="00:00")
  points(b$Lon[I],b$Lat[I],col=cm[i],pch=4,cex=2,lwd=1)
  a=unique(b$ID[I])
  for(j in 1:length(a)) 
  {
    c=b[b$ID==a[j],]
    lines(c$Lon,c$Lat,lwd=1,col=cm[i])
    points(c$Lon[c$Time=="00:00"],c$Lat[c$Time=="00:00"],lwd=1,col=cm[i],pch=4)
  }
  if(t==3) legend("topright",legend=names,lwd=1,col=cm,bg="white")
}
}
dev.off()

###Set 2 - Initial conditions

list=list(read.csv(paste(trackdir,"ECLfixes_default_31May.csv",sep="")),
          read.csv(paste(trackdir,"ECLfixes_default_fix.csv",sep="")),
          read.csv(paste(trackdir,"ECLfixes_default_2Jun.csv",sep="")))
names=c("31 May","1 Jun","2 Jun")
cm=c("red","blue","darkgreen")
keydates=c(20070609,20070619,20070627)

par(mar=c(4,4,2,2))
for(i in 1:3)
{
  b=list[[i]]
  d2=b$Date+(as.numeric(b$Time)-1)/4
  b$CV[b$Location==0]=NaN
  a=unique(b$ID)
  for(j in 1:length(a)) 
  {
    if(i==1 & j==1)
      plot(d2[b$ID==a[j]],b$CV[b$ID==a[j]],lwd=2,col=cm[i],type="l",
           xlim=c(20070601,20070630),ylim=c(0,4),
           xlab="Date",ylab="Curvature")
    else lines(d2[b$ID==a[j]],b$CV[b$ID==a[j]],lwd=2,col=cm[i])
  }
}
abline(v=keydates,col="gray")
legend("topleft",legend=names,lwd=2,col=cm)

par(mar=c(4,4,2,2))
for(i in 1:3)
{
  b=list[[i]]
  d2=b$Date+(as.numeric(b$Time)-1)/4
  b$MSLP[b$Location==0]=NaN
  a=unique(b$ID)
  for(j in 1:length(a)) 
  {
    if(i==1 & j==1)
      plot(d2[b$ID==a[j]],b$MSLP[b$ID==a[j]],lwd=2,col=cm[i],type="l",
           xlim=c(20070601,20070630),ylim=c(970,1030),
           xlab="Date",ylab="MSLP")
    else lines(d2[b$ID==a[j]],b$MSLP[b$ID==a[j]],lwd=2,col=cm[i])
  }
}
abline(v=keydates,col="gray")
legend("bottomleft",legend=names,lwd=2,col=cm)

pdf(file="WRFruns_Jun07_loc_initconds.pdf",width=8,height=3)
layout(cbind(1,2,3))
for(t in 1:length(keydates))
{
  par(mar=c(2,2,2,2))
  I=which(time==keydates[t]*100)
  contour(Useful$x,Useful$y,mask,drawlabels=F,xlim=c(145,180),ylim=c(-45,-10),main=paste(keydates[t],"00Z"))
  contour(lon,lat,slp[,,I],levels=seq(992,1050,4),add=T) #19980808 = 29 (1 + 7*4)
  for(i in 1:3)
  {
    b=list[[i]]
    I=which(b$Date==keydates[t] & b$Time=="00:00")
    points(b$Lon[I],b$Lat[I],col=cm[i],pch=4,cex=2,lwd=1)
    a=unique(b$ID[I])
    for(j in 1:length(a)) 
    {
      c=b[b$ID==a[j],]
      lines(c$Lon,c$Lat,lwd=1,col=cm[i])
      points(c$Lon[c$Time=="00:00"],c$Lat[c$Time=="00:00"],lwd=1,col=cm[i],pch=4)
    }
    if(t==3) legend("topright",legend=names,lwd=1,col=cm,bg="white")
  }
}
dev.off()

###All of 2007 tracks
setwd("~/Documents/ECLs")
library(RNetCDF)
library("R.matlab")
readMat('~/Documents/Data/Useful.mat')->Useful
mask<-t(Useful$mask)
mask[is.na(mask)]=0

trackdir="/home/nfs/z3478332/output/outputUM_wrf_2007/"
contour(Useful$x,Useful$y,mask,drawlabels=F,xlim=c(145,180),ylim=c(-45,-10))
control=read.csv(paste(trackdir,"ECLfixes_2007.csv",sep=""))
halftopo=read.csv(paste(trackdir,"ECLfixes_2007_halftopo.csv",sep=""))

lons=seq(145,180,2.5)
lats=seq(-45,-10,2.5)
cont<-half<-matrix(0,15,15)
for(i in 1:15)
  for(j in 1:15)
  {
    I=which(control$Lat>=lats[i]-1.25 & control$Lat<lats[i]+1.25 & 
              control$Lon>=lons[j]-1.25 & control$Lon<lons[j]+1.25)
    cont[i,j]=length(I)
    I=which(halftopo$Lat>=lats[i]-1.25 & halftopo$Lat<lats[i]+1.25 & 
              halftopo$Lon>=lons[j]-1.25 & halftopo$Lon<lons[j]+1.25)
    half[i,j]=length(I)
  }
library(fields)
library(s2dverification)

layout(cbind(c(1,3),c(2,3)),width=c(1,1),height=c(1,0.2))
image(lons,lats,cont,xlim=c(145,180),ylim=c(-45,-10),xlab="",ylab="",main="Default",
           breaks=seq(0,25,2.5),col=gray(seq(1,0.1,-0.1)),zlim=c(0,25))
contour(Useful$x,Useful$y,mask,drawlabels=F,add=T)
image(lons,lats,half,xlim=c(145,180),ylim=c(-45,-10),xlab="",ylab="",main="Halftopo",
           breaks=seq(0,25,2.5),col=gray(seq(1,0.1,-0.1)),zlim=c(0,25))
contour(Useful$x,Useful$y,mask,drawlabels=F,add=T)
ColorBar(brks=seq(0,25,2.5),cols=gray(seq(1,0.1,-0.1)), vert=FALSE, subsampleg = 2)



#####################Set 3 - SST changes


setwd("~/Documents/ECLs")
library(RNetCDF)
library("R.matlab")
readMat('~/Documents/Data/Useful.mat')->Useful
mask<-t(Useful$mask)
mask[is.na(mask)]=0

trackdir="/home/nfs/z3478332/output/outputUM_wrf_2007-06/"
wrfdir="/srv/ccrc/data36/z3478332/WRF/output/ERAI_R2_nudging_default/out/"
a=open.nc(paste(wrfdir,"WRF_d01_slp_regrid.nc",sep=""))
time=var.get.nc(a,"Times")
lat=var.get.nc(a,"lat0")
lon=var.get.nc(a,"lon0")
slp=var.get.nc(a,"slp0")

list=list(read.csv(paste(trackdir,"ECLfixes_default_SST-2.csv",sep="")),
          read.csv(paste(trackdir,"ECLfixes_default_SST-1.csv",sep="")),
          read.csv(paste(trackdir,"ECLfixes_default_fix.csv",sep="")),
          read.csv(paste(trackdir,"ECLfixes_default_SST+1.csv",sep="")),
          read.csv(paste(trackdir,"ECLfixes_default_SST+2.csv",sep="")))
names=c("SST -2K","SST -1K","Control","SST +1K","SST +2K")
cm=c("blue4","blue","black","red","red4")
keydates=c(20070609,20070619,20070627)

par(mar=c(4,4,2,2))
for(i in 1:5)
{
  b=list[[i]]
  d2=b$Date+(as.numeric(b$Time)-1)/4
  b$CV[b$Location==0]=NaN
  a=unique(b$ID)
  for(j in 1:length(a)) 
  {
    if(i==1 & j==1)
      plot(d2[b$ID==a[j]],b$CV[b$ID==a[j]],lwd=2,col=cm[i],type="l",
           xlim=c(20070601,20070630),ylim=c(0,4),
           xlab="Date",ylab="Curvature")
    else lines(d2[b$ID==a[j]],b$CV[b$ID==a[j]],lwd=2,col=cm[i])
  }
}
abline(v=keydates,col="gray")
legend("topleft",legend=names,lwd=2,col=cm)

par(mar=c(4,4,2,2))
for(i in 1:5)
{
  b=list[[i]]
  d2=b$Date+(as.numeric(b$Time)-1)/4
  b$MSLP[b$Location==0]=NaN
  a=unique(b$ID)
  for(j in 1:length(a)) 
  {
    if(i==1 & j==1)
      plot(d2[b$ID==a[j]],b$MSLP[b$ID==a[j]],lwd=2,col=cm[i],type="l",
           xlim=c(20070601,20070630),ylim=c(970,1030),
           xlab="Date",ylab="MSLP")
    else lines(d2[b$ID==a[j]],b$MSLP[b$ID==a[j]],lwd=2,col=cm[i])
  }
}
abline(v=keydates,col="gray")
legend("bottomleft",legend=names,lwd=2,col=cm)

pdf(file="WRFruns_Jun07_loc_SST.pdf",width=8,height=3)
layout(cbind(1,2,3))
for(t in 1:length(keydates))
{
  par(mar=c(2,2,2,2))
  I=which(time==keydates[t]*100)
  contour(Useful$x,Useful$y,mask,drawlabels=F,xlim=c(145,180),ylim=c(-45,-10),main=paste(keydates[t],"00Z"))
  contour(lon,lat,slp[,,I],levels=seq(992,1050,4),add=T) #19980808 = 29 (1 + 7*4)
  for(i in 1:5)
  {
    b=list[[i]]
    I=which(b$Date==keydates[t] & b$Time=="00:00")
    points(b$Lon[I],b$Lat[I],col=cm[i],pch=4,cex=2,lwd=1)
    a=unique(b$ID[I])
    for(j in 1:length(a)) 
    {
      c=b[b$ID==a[j],]
      lines(c$Lon,c$Lat,lwd=1,col=cm[i])
      points(c$Lon[c$Time=="00:00"],c$Lat[c$Time=="00:00"],lwd=1,col=cm[i],pch=4)
    }
    if(t==3) legend("topright",legend=names,lwd=1,col=cm,bg="white")
  }
}
dev.off()



#####################Set 4 - WRF Versions


setwd("~/Documents/ECLs")
library(RNetCDF)
library("R.matlab")
readMat('~/Documents/Data/Useful.mat')->Useful
mask<-t(Useful$mask)
mask[is.na(mask)]=0

trackdir="/home/nfs/z3478332/output/outputUM_wrf_2007-06/"
wrfdir="/srv/ccrc/data36/z3478332/WRF/output/ERAI_R2_nudging_default/out/"
a=open.nc(paste(wrfdir,"WRF_d01_slp_regrid.nc",sep=""))
time=var.get.nc(a,"Times")
lat=var.get.nc(a,"lat0")
lon=var.get.nc(a,"lon0")
slp=var.get.nc(a,"slp0")

list=list(read.csv(paste(trackdir,"ECLfixes_R1_nudging_default.csv",sep="")),
          read.csv(paste(trackdir,"ECLfixes_default_fix.csv",sep="")),
          read.csv(paste(trackdir,"ECLfixes_R3_nudging_default.csv",sep="")),
          read.csv(paste(trackdir,"ECLfixes_R1_nudging_default_50.csv",sep="")),
          read.csv(paste(trackdir,"ECLfixes_default_50.csv",sep="")),
          read.csv(paste(trackdir,"ECLfixes_R3_nudging_default_50.csv",sep="")))
names=c("R1","R2","R3")
cm=c("red","blue","darkgreen")
keydates=c(20070609,20070619,20070627)

par(mar=c(4,4,2,2))
for(i in 1:3)
{
  b=list[[i]]
  d2=b$Date+(as.numeric(b$Time)-1)/4
  b$CV[b$Location==0]=NaN
  a=unique(b$ID)
  for(j in 1:length(a)) 
  {
    if(i==1 & j==1)
      plot(d2[b$ID==a[j]],b$CV[b$ID==a[j]],lwd=2,col=cm[i],type="l",
           xlim=c(20070601,20070630),ylim=c(0,4),
           xlab="Date",ylab="Curvature")
    else lines(d2[b$ID==a[j]],b$CV[b$ID==a[j]],lwd=2,col=cm[i])
  }
}
abline(v=keydates,col="gray")
legend("topleft",legend=names,lwd=2,col=cm)

par(mar=c(4,4,2,2))
for(i in 1:3)
{
  b=list[[i]]
  d2=b$Date+(as.numeric(b$Time)-1)/4
  b$MSLP[b$Location==0]=NaN
  a=unique(b$ID)
  for(j in 1:length(a)) 
  {
    if(i==1 & j==1)
      plot(d2[b$ID==a[j]],b$MSLP[b$ID==a[j]],lwd=2,col=cm[i],type="l",
           xlim=c(20070601,20070630),ylim=c(970,1030),
           xlab="Date",ylab="MSLP")
    else lines(d2[b$ID==a[j]],b$MSLP[b$ID==a[j]],lwd=2,col=cm[i])
  }
}
abline(v=keydates,col="gray")
legend("bottomleft",legend=names,lwd=2,col=cm)

pdf(file="WRFruns_Jun07_loc_WRFVersion.pdf",width=8,height=3)
layout(cbind(1,2,3))
for(t in 1:length(keydates))
{
  par(mar=c(2,2,2,2))
  I=which(time==keydates[t]*100)
  contour(Useful$x,Useful$y,mask,drawlabels=F,xlim=c(145,180),ylim=c(-45,-10),main=paste(keydates[t],"00Z"))
  contour(lon,lat,slp[,,I],levels=seq(992,1050,4),add=T) #19980808 = 29 (1 + 7*4)
  for(i in 1:3)
  {
    b=list[[i]]
    I=which(b$Date==keydates[t] & b$Time=="00:00")
    points(b$Lon[I],b$Lat[I],col=cm[i],pch=4,cex=2,lwd=1)
    a=unique(b$ID[I])
    for(j in 1:length(a)) 
    {
      c=b[b$ID==a[j],]
      lines(c$Lon,c$Lat,lwd=1,col=cm[i])
      points(c$Lon[c$Time=="00:00"],c$Lat[c$Time=="00:00"],lwd=1,col=cm[i],pch=4)
    }
    if(t==3) legend("topright",legend=names,lwd=1,col=cm,bg="white")
  }
}
dev.off()

par(mar=c(4,4,2,2))
for(i in 1:3)
{
  b=list[[i]]
  d2=b$Date+(as.numeric(b$Time)-1)/4
  b$CV[b$Location==0]=NaN
  a=unique(b$ID)
  for(j in 1:length(a)) 
  {
    if(i==1 & j==1)
      plot(d2[b$ID==a[j]],b$CV[b$ID==a[j]],lwd=2,col=cm[i],type="l",
           xlim=c(20070601,20070630),ylim=c(0,4),
           xlab="Date",ylab="Curvature")
    else lines(d2[b$ID==a[j]],b$CV[b$ID==a[j]],lwd=2,col=cm[i])
  }
  b=list[[i+3]]
  d2=b$Date+(as.numeric(b$Time)-1)/4
  b$CV[b$Location==0]=NaN
  a=unique(b$ID)
  for(j in 1:length(a)) lines(d2[b$ID==a[j]],b$CV[b$ID==a[j]],lty=2,lwd=2,col=cm[i])
}
abline(v=keydates,col="gray")
legend("topleft",legend=names,lwd=2,col=cm)

par(mar=c(4,4,2,2))
for(i in 1:3)
{
  b=list[[i]]
  d2=b$Date+(as.numeric(b$Time)-1)/4
  b$MSLP[b$Location==0]=NaN
  a=unique(b$ID)
  for(j in 1:length(a)) 
  {
    if(i==1 & j==1)
      plot(d2[b$ID==a[j]],b$MSLP[b$ID==a[j]],lwd=2,col=cm[i],type="l",
           xlim=c(20070601,20070630),ylim=c(970,1030),
           xlab="Date",ylab="MSLP")
    else lines(d2[b$ID==a[j]],b$MSLP[b$ID==a[j]],lwd=2,col=cm[i])
  }
  b=list[[i+3]]
  d2=b$Date+(as.numeric(b$Time)-1)/4
  b$MSLP[b$Location==0]=NaN
  a=unique(b$ID)
  for(j in 1:length(a))  lines(d2[b$ID==a[j]],b$MSLP[b$ID==a[j]],lty=2,lwd=2,col=cm[i])
}
abline(v=keydates,col="gray")
legend("bottomleft",legend=names,lwd=2,col=cm)

pdf(file="WRFruns_Jun07_loc_WRFVersion_halftopo.pdf",width=8,height=3)
layout(cbind(1,2,3))
for(t in 1:length(keydates))
{
  par(mar=c(2,2,2,2))
  I=which(time==keydates[t]*100)
  contour(Useful$x,Useful$y,mask,drawlabels=F,xlim=c(145,180),ylim=c(-45,-10),main=paste(keydates[t],"00Z"))
  contour(lon,lat,slp[,,I],levels=seq(992,1050,4),add=T) #19980808 = 29 (1 + 7*4)
  for(i in 1:3)
  {
    b=list[[i]]
    I=which(b$Date==keydates[t] & b$Time=="00:00")
    points(b$Lon[I],b$Lat[I],col=cm[i],pch=4,cex=2,lwd=1)
    a=unique(b$ID[I])
    for(j in 1:length(a)) 
    {
      c=b[b$ID==a[j],]
      lines(c$Lon,c$Lat,lwd=1,col=cm[i])
      points(c$Lon[c$Time=="00:00"],c$Lat[c$Time=="00:00"],lwd=1,col=cm[i],pch=4)
    }
    b=list[[i+3]]
    I=which(b$Date==keydates[t] & b$Time=="00:00")
    a=unique(b$ID[I])
    for(j in 1:length(a)) 
    {
      c=b[b$ID==a[j],]
      lines(c$Lon,c$Lat,lwd=1,col=cm[i],lty=2)
    }
    if(t==3) legend("topright",legend=names,lwd=1,col=cm,bg="white")
  }
}
dev.off()
