###Looking at a u-wind based index?
rm(list=ls())
setwd('~/Documents/Data/JE_WRF')
library("RNetCDF")

source('~/Documents/R/filled.contour3.R')
source('~/Documents/R/filled.legend.R')##Plotty stuff   
source('~/Documents/R/color.palette.R')
library("R.matlab")
readMat('~/Documents/Data/Useful.mat')->Useful
mask<-t(Useful$mask)
mask[is.na(mask)]=0
pal <- color.palette(c("red","yellow","white","cyan","blue"), c(20,5,5,20))
cm=pal(20)
cm[9:12]="white"

f1=open.nc('~/Documents/Data/NCEP Uwind/uwnd.mon.mean.nc')
time=var.get.nc(f1,'time')
lat=var.get.nc(f1,'lat')
lon=var.get.nc(f1,'lon')
time=as.Date((time/24-2),origin="0001-01-01")
time=cbind(as.numeric(format(time,"%Y")),as.numeric(format(time,"%m")))
uwnd=var.get.nc(f1,'uwnd',c(1,1,3,1),c(144,73,1,length(time[,1])),unpack=T)
I=which(time[,1]>=1979 & time[,1]<=2009)
uwnd2=uwnd[,,I]
Uind1=apply(uwnd2[61,47:51,],2,mean)
Uind2=apply(uwnd2[59:61,47:51,],3,mean)
Uwind=array(0,c(31,12,2))
for(i in 1:length(years))
  for(j in 1:length(months))
  {
    I=which(time2[,2]==years[i] & time2[,3]==j)
    Uwind[i,j,1]=Uind1[I]
    Uwind[i,j,2]=Uind2[I]
  }

source('~/Documents/R/WRF_corrs.R')
seasons=c("ann","warm","cool","mam","jja","son","djf")
snum=cbind(c(1,11,5,3,6,9,12),c(12,4,10,5,8,11,2))
for(i in 1:7)
{
  gen.corrA(snum[i,],seasons[i],Uwind[,,1],'Uwind1')
  gen.corrA(snum[i,],seasons[i],Uwind[,,2],'Uwind2')
}

save(Uind1,Uind2,Uwind,file="Uwind_NCEP.RData")

##Okay, so now we can load in the WRF U-winds!
##Re-do with better WRF
rm(list=ls())
setwd('~/Documents/Data/JE_WRF')
library("RNetCDF")

f1=open.nc('~/Documents/Data/WRF_meanwind.nc')
U=var.get.nc(f1,'U')
lat=var.get.nc(f1,'lat')
Uave=apply(U[111:114,49:70,],3,mean,na.rm=T)
load('Uwind.RData')

comp=matrix(0,12,6)
for(i in 1:12)
{
  I=which(Uwind[,2]==i)
  comp[i,1]=mean(Uwind[I,3])
  comp[i,2]=mean(Uwind[I,4])
  comp[i,3]=mean(Uave[I])
  comp[i,4]=cor(Uwind[I,3],Uwind[I,4])
  comp[i,5]=cor(Uwind[I,3],Uave[I])
  comp[i,6]=cor(Uwind[I,4],Uave[I])
}
mnames=c("Jan","Feb","Mar","Apr","May","Jun","Jul","Aug","Sep","Oct","Nov","Dec")
plot(seq(1:12),comp[,1],ty="l",lwd=2,col="blue",xlab="",ylab="",axes=F,ylim=c(-3,4))
axis(1, at = seq(1:12), labels = c(mnames), cex.axis = 0.75) 
axis(2, at = seq(-3,4))
title(main="Average Uwind (m/s)",cex.main=1)
lines(seq(1:12),comp[,3],ty="l",lwd=2,col="red")
abline(h=0,col="gray",lty=2)
legend("topleft",legend=c("NCEP","WRF"),col=c("blue","red"),lwd=2,bty="n")

##
read.table('GDI.csv',header=T,sep=",")->gdi

##Now, do Uwind/rain corrs for all seasons for WRF
years=seq(1979,2009)
U=matrix(0,31,12)
for(i in 1:31)
  for(j in 1:12)
  {
    I=which(Uwind[,1]==years[i] & Uwind[,2]==j)
    U[i,j]=Uave[I]    
  }

##
rm(list=ls())
setwd('~/Documents/Data/JE_WRF')
library("RNetCDF")
source('~/Documents/R/WRF_corrs.R')
seasons=c("ann","warm","cool","mam","jja","son","djf")
snum=cbind(c(1,11,5,3,6,9,12),c(12,4,10,5,8,11,2))

f1=open.nc('~/Documents/Data/WRF_meanwind.nc')
U=var.get.nc(f1,'U')
lat=var.get.nc(f1,'lat')
Uave=apply(U[111:114,37:59,],3,mean,na.rm=T) ##40-30S
U=matrix(0,31,12)
for(i in 1:31)
  for(j in 1:12)
  {
    I=which(time2[,2]==years[i] & time2[,3]==j)
    U[i,j]=Uave[I]    
  }
for(i in 1:7) gen.corr2(snum[i,],seasons[i],U,'UwindW_S')

f1=open.nc('~/Documents/Data/NCEP Uwind/uwnd.mon.mean.nc')
time=var.get.nc(f1,'time')
time=as.Date((time/24-2),origin="0001-01-01")
time=cbind(as.numeric(format(time,"%Y")),as.numeric(format(time,"%m")))
uwnd=var.get.nc(f1,'uwnd',c(1,1,3,1),c(144,73,1,length(time[,1])),unpack=T)
I=which(time[,1]>=1979 & time[,1]<=2009)
uwnd2=uwnd[,,I]
Uind1=apply(uwnd2[61,49:53,],2,mean)
U=matrix(0,31,12)
for(i in 1:31)
  for(j in 1:12)
  {
    I=which(time2[,2]==years[i] & time2[,3]==j)
    U[i,j]=Uind1[I]    
  }
for(i in 1:7) gen.corr2(snum[i,],seasons[i],U,'UwindN_S')

###Daily Uwind
rm(list=ls())
setwd('~/Documents/Data/JE_WRF')
library("RNetCDF")
read.csv('WRF_Uwind.csv',sep=";")->WRFdaily

##Want lats 11-15, lon 21
years=seq(1979,2009)
for(i in 1:length(years))
{
  fin=paste('~/Documents/Data/NCEP Uwind/uwnd.850.',years[i],'.nc',sep="")
  f1=open.nc(fin)
  time=var.get.nc(f1,'time')
  uwnd=var.get.nc(f1,'uwnd',c(21,11,1,1),c(1,5,1,length(time)),unpack=T)
  hh=(time %% 24)
  if(i==1) NCEPdaily=apply(uwnd[,hh==0],2,mean,na.rm=T) else NCEPdaily=c(NCEPdaily,apply(uwnd[,hh==0],2,mean,na.rm=T) )
}

##Prop E/W
Uwind=matrix(0,31*12,4)
n=1
for(i in 1:length(years))
  for(j in 1:12)
  {
    I=which(WRFdaily[,2]==years[i] & WRFdaily[,3]==j)
    Uwind[n,1]=years[i]
    Uwind[n,2]=j  
    J=which(WRFdaily[I,5]>0)
    Uwind[n,3]=length(J)/length(I)
    J=which(NCEPdaily[I]>0)
    Uwind[n,4]=length(J)/length(I)
    n=n+1
  }
comp=matrix(0,12,3)
for(i in 1:12)
{
  I=which(Uwind[,2]==i)
  comp[i,1]=mean(Uwind[I,3])
  comp[i,2]=mean(Uwind[I,4])
  comp[i,3]=cor(Uwind[I,3],Uwind[I,4])
}

###############
##Prop rain on E/W days
##Using both AWAP & WRF
source('~/Documents/R/WRF_corrs.R')
load('Rain_daily.RData')
seasons=c("ann","warm","cool","mam","jja","son","djf")
snum=cbind(c(1,11,5,3,6,9,12),c(12,4,10,5,8,11,2))
for(i in 1:7) gen.propW(snum[i,],seasons[i],rain2,-WRFdaily[,5],'Uwrf')
for(i in 1:7) gen.propW(snum[i,],seasons[i],rain2,-NCEPdaily,'Uncep')
for(i in 1:7) gen.propA(snum[i,],seasons[i],-NCEPdaily,'Uncep')