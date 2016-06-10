
##Currently edited to produce figs for WRF10 not WRF50, but quick to change refs.
rm(list=ls())
setwd("~/Documents/Data/CMIP5/CMIP3")
library(RNetCDF)

name=c("MPI_ECHAM5","CSIRO_Mk3.0","MIROC3.2","CGCM3.1")
Ufile=c("/srv/ccrc/data34/z3478332/CMIP3/uw850/hus_850hPa_picntrl_mpi-echam5_2.5deg.nc",
        "/srv/ccrc/data34/z3478332/CMIP3/uw850/hus_850hPa_picntrl_csiromk3_2.5deg.nc",
        "/srv/ccrc/data34/z3478332/CMIP3/uw850/hus_850hPa_picntrl_miroc3.2_2.5deg.nc",
        "/srv/ccrc/data34/z3478332/CMIP3/uw850/hus_850hPa_picntrl_cgcm3.1_2.5deg.nc")
Wfile=c("/srv/ccrc/data34/z3478332/CMIP3/WRF_echam5_monthly.nc",
        "/srv/ccrc/data34/z3478332/CMIP3/WRF_csiromk3_monthly.nc",
        "/srv/ccrc/data34/z3478332/CMIP3/WRF_miroc3.2_monthly.nc",
        "/srv/ccrc/data34/z3478332/CMIP3/WRF_cccma_monthly.nc")
Wfile10=c("/srv/ccrc/data34/z3478332/CMIP3/WRF10_echam5_monthly.nc",
        "/srv/ccrc/data34/z3478332/CMIP3/WRF10_csiromk3_monthly.nc",
        "/srv/ccrc/data34/z3478332/CMIP3/WRF10_miroc3.2_monthly.nc",
        "/srv/ccrc/data34/z3478332/CMIP3/WRF10_cccma_monthly.nc")

##Need to have NCEP for comparison! 
a=23
b=47
f1=open.nc('~/Documents/Data/NCEP Uwind/uwnd.mon.mean.nc')
timeN=var.get.nc(f1,'time')
timeN=as.Date((timeN/24-2),origin="0001-01-01")

timeN=cbind(as.numeric(format(timeN,"%Y")),as.numeric(format(timeN,"%m")))
uwnd=var.get.nc(f1,'uwnd',c(1,1,3,1),c(144,73,1,length(timeN[,1])),unpack=T)
NCEP=apply(var.get.nc(f1,'uwnd',c(61,b,3,1),c(1,5,1,length(timeN[,1])),unpack=T),2,mean,na.rm=T)

for(n in 1:4)
{
print(name[n])
  fig1=paste("Uwind_",name[n],"W10.tiff",sep="")
fig2a=paste("Rain_corrUwind_sig_",name[n],"W10_cool.tiff",sep="")
fig3a=paste("Rain_corrUwind_sig_",name[n],"W10_warm.tiff",sep="")
fig2=paste("Rain_corrUwind_",name[n],"W10_cool.tiff",sep="")
fig3=paste("Rain_corrUwind_",name[n],"W10_warm.tiff",sep="")

f1=open.nc(Wfile10[n])
rain=var.get.nc(f1,"prcp")
uwind=var.get.nc(f1,"uw850")
lat=var.get.nc(f1,"lat")
lon=var.get.nc(f1,"lon")
time=var.get.nc(f1,"time")
year=floor(time/100)
month=time%%100
time=cbind(year,month)
Uave=apply(uwind[61,a:(a+4),],2,mean)

##Actually, let's now include the CMIP3 version
f1=open.nc(Ufile[n])
u=var.get.nc(f1,"ua")
Uave2=apply(u[61,a:(a+4),],2,mean)
timeC=var.get.nc(f1,"time")
ori=att.get.nc(f1,"time","units")
ori=unlist(strsplit(ori, split=" "))[3]
time2=as.Date(timeC,origin=ori)
timeC=data.frame(Year=as.numeric(format(time2,"%Y")),Month=as.numeric(format(time2,"%m")))

if(n==2) {
#Fix for CSIRO
yy=seq(1931,2030)
mm=seq(1,12)
k=1
for(i in 1:100)
  for(j in 1:12)
  {
    timeC[k,1]=yy[i]
    timeC[k,2]=mm[j]
    k=k+1
  }
} 

if(n==3) I=which(timeC[,1]>=80 & timeC[,1]<=100) else I=which(timeC[,1]>=1990 & timeC[,1]<=2009)
Uave2=Uave2[I]
timeC=timeC[I,]

##Edited version w/ WRF10 & WRF50
f1=open.nc(Wfile[n])
uwind50=var.get.nc(f1,"uw850")
Uave50=apply(uwind50[61,a:(a+4),],2,mean)

Ucomp=matrix(0,12,5)
for(i in 1:12)
{
  Ucomp[i,1]=i
  I=which(timeN[,1]<=2009 & timeN[,1]>=1990 & timeN[,2]==i)
  Ucomp[i,2]=mean(NCEP[I])
  I=which(time[,2]==i)
  Ucomp[i,4]=mean(Uave50[I])
  Ucomp[i,5]=mean(Uave[I])
  I=which(timeC[,2]==i)
  Ucomp[i,3]=mean(Uave2[I])
}

mnames=c("Jan","Feb","Mar","Apr","May","Jun","Jul","Aug","Sep","Oct","Nov","Dec")
yl=c(floor(min(Ucomp[,2:3])),ceiling(max(Ucomp[,2:3])))
tiff(file=fig1, height=500, width=800)
plot(seq(1:12),Ucomp[,2],ty="l",lwd=2,col="blue",xlab="",ylab="",axes=F,ylim=yl)
axis(1, at = seq(1:12), labels = c(mnames),cex.axis=1.5) 
axis(2, at = seq(yl[1],yl[2]),cex.axis=1.5)
title(main="Average Uwind (m/s)",cex.main=2)
lines(seq(1:12),Ucomp[,3],ty="l",lwd=2,col="red")
lines(seq(1:12),Ucomp[,4],ty="l",lwd=2,col="green")
lines(seq(1:12),Ucomp[,5],ty="l",lwd=2,col="purple")
abline(h=0,col="gray",lty=2)
legend("topleft",legend=c("NCEP",name[n],"WRF 50km","WRF 10km"),col=c("blue","red","green","purple"),lwd=2,bty="n",cex=1.5)
dev.off()

#Finally, let's do the two correlation plots
years=seq(min(time[,1]),max(time[,1]))
years=cbind(seq(min(time[,1]),max(time[,1])),matrix(0,(length(years)),2))
coolR<-warmR<-array(0,c(144,73,length(years[,1])))
for(i in 1:length(years[,1]))
{
  I=which(time[,1]==years[i,1] & time[,2]>=5 & time[,2]<=10)
  years[i,2]=mean(Uave[I])
  coolR[,,i]=apply(rain[,,I],c(1,2),sum)
  I=which((time[,1]==years[i,1] & time[,2]<=4) | (time[,1]==years[i,1]-1 & time[,2]>=11))
  years[i,3]=mean(Uave[I])
  warmR[,,i]=apply(rain[,,I],c(1,2),sum)
}

coolC<-warmC<-matrix(0,144,73)
for(i in 1:144)
  for(j in 1:73)
  {
    coolC[i,j]=cor(coolR[i,j,],years[,2],use="pairwise.complete.obs")
    warmC[i,j]=cor(warmR[i,j,],years[,3],use="pairwise.complete.obs")
  }

i=which(lat==-35)
j=which(lon==145)
print(coolC[j,i])
i=which(lat==-30)
j=which(lon==152.5)
print(coolC[j,i])

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

coolC[coolC>=0.69]=0.69
coolC[coolC<=(-0.69)]=-0.69
warmC[warmC>=0.69]=0.69
warmC[warmC<=(-0.69)]=-0.69

tiff(file=fig2, height=500, width=800)
image.plot(lon,lat,coolC,xlab="",breaks=seq(-0.7,0.7,0.1),ylab="",col=cm,zlim=c(-0.7,0.7),xlim=c(100,180),ylim=c(-50,-10))
contour(Useful$x,Useful$y,mask,add=T,drawlabels=F)
dev.off()
tiff(file=fig3, height=500, width=800)
image.plot(lon,lat,warmC,xlab="",breaks=seq(-0.7,0.7,0.1),ylab="",col=cm,zlim=c(-0.7,0.7),xlim=c(100,180),ylim=c(-50,-10))
contour(Useful$x,Useful$y,mask,add=T,drawlabels=F)
dev.off()

coolC[abs(coolC)<=0.27]=NaN
warmC[abs(warmC)<=0.27]=NaN

tiff(file=fig2a, height=500, width=800)
image.plot(lon,lat,coolC,xlab="",breaks=seq(-0.7,0.7,0.1),ylab="",col=cm,zlim=c(-0.7,0.7),xlim=c(100,180),ylim=c(-50,-10))
contour(Useful$x,Useful$y,mask,add=T,drawlabels=F)
dev.off()
tiff(file=fig3a, height=500, width=800)
image.plot(lon,lat,warmC,xlab="",breaks=seq(-0.7,0.7,0.1),ylab="",col=cm,zlim=c(-0.7,0.7),xlim=c(100,180),ylim=c(-50,-10))
contour(Useful$x,Useful$y,mask,add=T,drawlabels=F)
dev.off()
}


##Comparing the four versions for MEDIANS etc.

WRF_compW<-WRF_compC<-WRF10_compW<-WRF10_compC<-array(NaN,c(144,73,4))

for(n in 1:4)
{
  f1=open.nc(Wfile10[n])
  rain=var.get.nc(f1,"prcp")
  uwind=var.get.nc(f1,"uw850")
  lat=var.get.nc(f1,"lat")
  lon=var.get.nc(f1,"lon")
  time=var.get.nc(f1,"time")
  year=floor(time/100)
  month=time%%100
  time=cbind(year,month)
  Uave=apply(uwind[61,a:(a+4),],2,mean)
  
  years=seq(min(time[,1]),max(time[,1]))
  years=cbind(seq(min(time[,1]),max(time[,1])),matrix(0,(length(years)),2))
  coolR<-warmR<-array(0,c(144,73,length(years[,1])))
  for(i in 1:length(years[,1]))
  {
    I=which(time[,1]==years[i,1] & time[,2]>=5 & time[,2]<=10)
    years[i,2]=mean(Uave[I])
    coolR[,,i]=apply(rain[,,I],c(1,2),sum)
    I=which((time[,1]==years[i,1] & time[,2]<=4) | (time[,1]==years[i,1]-1 & time[,2]>=11))
    years[i,3]=mean(Uave[I])
    warmR[,,i]=apply(rain[,,I],c(1,2),sum)
  }
  
  coolC<-warmC<-matrix(0,144,73)
  for(i in 1:144)
    for(j in 1:73)
    {
      WRF10_compC[i,j,n]=cor(coolR[i,j,],years[,2],use="pairwise.complete.obs")
      WRF10_compW[i,j,n]=cor(warmR[i,j,],years[,3],use="pairwise.complete.obs")
    }
  
  f1=open.nc(Wfile[n])
  rain=var.get.nc(f1,"prcp")
  uwind=var.get.nc(f1,"uw850")
  lat=var.get.nc(f1,"lat")
  lon=var.get.nc(f1,"lon")
  time=var.get.nc(f1,"time")
  year=floor(time/100)
  month=time%%100
  time=cbind(year,month)
  Uave=apply(uwind[61,a:(a+4),],2,mean)
  
  years=seq(min(time[,1]),max(time[,1]))
  years=cbind(seq(min(time[,1]),max(time[,1])),matrix(0,(length(years)),2))
  coolR<-warmR<-array(0,c(144,73,length(years[,1])))
  for(i in 1:length(years[,1]))
  {
    I=which(time[,1]==years[i,1] & time[,2]>=5 & time[,2]<=10)
    years[i,2]=mean(Uave[I])
    coolR[,,i]=apply(rain[,,I],c(1,2),sum)
    I=which((time[,1]==years[i,1] & time[,2]<=4) | (time[,1]==years[i,1]-1 & time[,2]>=11))
    years[i,3]=mean(Uave[I])
    warmR[,,i]=apply(rain[,,I],c(1,2),sum)
  }
  
  coolC<-warmC<-matrix(0,144,73)
  for(i in 1:144)
    for(j in 1:73)
    {
      WRF_compC[i,j,n]=cor(coolR[i,j,],years[,2],use="pairwise.complete.obs")
      WRF_compW[i,j,n]=cor(warmR[i,j,],years[,3],use="pairwise.complete.obs")
    }
}

medW=apply(WRF10_compW,c(1,2),median,na.rm=T)
medW[medW>=0.69]=0.69
medW[medW<=(-0.69)]=-0.69
#medW[abs(medW)<=0.27]=NaN
medC=apply(WRF10_compC,c(1,2),median,na.rm=T)
medC[medC>=0.69]=0.69
medC[medC<=(-0.69)]=-0.69
#medC[abs(medC)<=0.27]=NaN

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
tiff(file=paste("~/Documents/Data/CMIP5/Rain_corrUwind_WRF10R2_median_cool_small.tiff",sep=""), height=500, width=600,pointsize=20)
image.plot(lon,lat,medC,xlab="",breaks=seq(-0.7,0.7,0.1),ylab="",col=cm,zlim=c(-0.7,0.7),xlim=c(110,160),ylim=c(-45,-10))
contour(Useful$x,Useful$y,mask,add=T,drawlabels=F)
dev.off()
tiff(file=paste("~/Documents/Data/CMIP5/Rain_corrUwind_WRF10R2_median_warm_small.tiff",sep=""), height=500, width=600,pointsize=20)
image.plot(lon,lat,medW,xlab="",breaks=seq(-0.7,0.7,0.1),ylab="",col=cm,zlim=c(-0.7,0.7),xlim=c(110,160),ylim=c(-45,-10))
contour(Useful$x,Useful$y,mask,add=T,drawlabels=F)
dev.off()

medW[abs(medW)<=0.27]=NaN
medC[abs(medC)<=0.27]=NaN
tiff(file=paste("~/Documents/Data/CMIP5/Rain_corrUwind_sig_WRF10R2_median_cool_small.tiff",sep=""), height=500, width=600,pointsize=20)
image.plot(lon,lat,medC,xlab="",breaks=seq(-0.7,0.7,0.1),ylab="",col=cm,zlim=c(-0.7,0.7),xlim=c(110,160),ylim=c(-45,-10))
contour(Useful$x,Useful$y,mask,add=T,drawlabels=F)
dev.off()
tiff(file=paste("~/Documents/Data/CMIP5/Rain_corrUwind_sig_WRF10R2_median_warm_small.tiff",sep=""), height=500, width=600,pointsize=20)
image.plot(lon,lat,medW,xlab="",breaks=seq(-0.7,0.7,0.1),ylab="",col=cm,zlim=c(-0.7,0.7),xlim=c(110,160),ylim=c(-45,-10))
contour(Useful$x,Useful$y,mask,add=T,drawlabels=F)
dev.off()











#############
##Wind comps, for reference

library(RNetCDF)
f1=open.nc("/srv/ccrc/data34/z3478332/CMIP3/WRF_R2_csiromk3_uwnd.nc")
uwind=var.get.nc(f1,"uwnd250")
a=23
Uave1=apply(uwind[61,a:(a+4),],2,mean)

uwind=var.get.nc(f1,"uwnd50")
lat=var.get.nc(f1,"lat50")
lon=var.get.nc(f1,"lon50")
time=var.get.nc(f1,"time")
I=which(lat>=-36.25 & lat<=-23.75)
J=which(lon>=148.75 & lon<=151.25)
Uave2=apply(uwind[J,I,],3,mean)

date=cbind(floor(time/100),time%%100)

Ucomp2=matrix(0,13,5)
Ucomp2[,1]=seq(0,12)
Ucomp2[1,2]=mean(Uave1)
Ucomp2[1,3]=mean(Uave2)
Ucomp2[1,5]=cor(Uave1,Uave2)
for(i in 1:12)
{
  I=which(date[,2]==i)
  Ucomp2[i+1,2]=mean(Uave1[I])
  Ucomp2[i+1,3]=mean(Uave2[I])
  Ucomp2[i+1,5]=cor(Uave1[I],Uave2[I])
}
Ucomp2[,4]=Ucomp2[,2]-Ucomp2[,3]

plot(Ucomp2[2:13,1],Ucomp2[2:13,2],type="l",col="red",ylim=c(-5,10))
lines(Ucomp2[2:13,1],Ucomp2[2:13,3],type="l",col="blue")
legend("topright",c("2.5 degree","0.5 degree"),col=c("red","blue"),lwd=c(1,1))
