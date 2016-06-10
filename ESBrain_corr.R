###Looking at a u-wind based index?
rm(list=ls())
setwd('~/Documents/Data/JE_WRF')
library("RNetCDF")
source('~/Documents/R/WRF_corrs.R')

##Which has the strongest correlation with ESB rain?
readMat('~/Documents/Data/mask_escci2.mat')->escci
escci<-t(escci$mask)
ESBrain=rep(0,372)
for(i in 1:372) ESBrain[i]=mean(AWAP[,,i]*escci,na.rm=T)
seasons=c("ann","warm","cool","mam","jja","son","djf")
smon=cbind(c(1,11,5,3,6,9,12),c(12,4,10,5,8,11,2))

##Using NCEP
f1=open.nc('~/Documents/Data/uwnd.mon.mean.nc')
time=var.get.nc(f1,'time')
lat=var.get.nc(f1,'lat')
lon=var.get.nc(f1,'lon')
time=as.Date((time/24-2),origin="0001-01-01")
time=cbind(as.numeric(format(time,"%Y")),as.numeric(format(time,"%m")))
uwnd=var.get.nc(f1,'uwnd',c(1,1,3,1),c(144,73,1,length(time[,1])),unpack=T)
I=which(time[,1]>=1979 & time[,1]<=2009)
uwnd2=uwnd[,,I]
for(i in 1:7) esb.corr(smon[i,],seasons[i],uwnd2,'UwindN','AWAP',ESBrain,lat,lon)

##Now for MSLP-based
f1=open.nc('~/Documents/Data/ncep.mslp.nc')
time=var.get.nc(f1,'time')
lat=var.get.nc(f1,'lat')
lon=var.get.nc(f1,'lon')
time=as.Date((time/24-2),origin="0001-01-01")
time=cbind(as.numeric(format(time,"%Y")),as.numeric(format(time,"%m")))
mslp=var.get.nc(f1,'slp',c(1,1,1),c(144,73,length(time[,1])),unpack=T)
I=which(time[,1]>=1979 & time[,1]<=2009)
slp2=mslp[,,I]
for(i in 1:7) esb.corr(smon[i,],seasons[i],slp2,'MSLP','AWAP',ESBrain,lat,lon)

##Next, need a WRF equivalent of ESB rain
##Best is probably to convert WRF to a 691x886 matrix then do as for AWAP

ESBW=rep(0,372)
for(i in 1:372)
{
  a=as.vector(monthR[,,i])
  a[which(is.na(a) | is.infinite(a))]=0
  b=interp(lont,latt,a,Useful$x,Useful$y)
  ESBW[i]=mean(b$z*escci,na.rm=T)
}

save(ESBW,ESBrain,file="ESBrain.RData")

f1=open.nc('~/Documents/Data/WRF_meanwind.nc')
U=var.get.nc(f1,'U')
lat=var.get.nc(f1,'lat')
lon=var.get.nc(f1,'lon')
for(i in 1:7) esb.corr(smon[i,],seasons[i],U,'UwindW','WRF',ESBW,lat,lon)

years=seq(1979,2009)
WRF<-array(0,dim=c(215,144,372))
n=1
for(i in 1:length(years))
  for(m in 1:12)
  {
    if(m<10) fname=paste('/srv/ccrc/data23/z3444417/studies/Data/WRF50/nnrp/R2/psl/raw/monthly/WRF50_R2_slp_',years[i],'-0',m,'.nc',sep="")
    else fname=paste('/srv/ccrc/data23/z3444417/studies/Data/WRF50/nnrp/R2/psl/raw/monthly/WRF50_R2_slp_',years[i],'-',m,'.nc',sep="")
    f1=open.nc(fname)
    nt=dim.inq.nc(f1,'Time')
    nt<-nt$length
    slp=var.get.nc(f1,'SLP',c(1,1,1),c(215,144,nt))
    close.nc(f1)
    WRF[,,n]=apply(slp[,,],c(1,2),mean,na.rm=T)
    n=n+1
  }
save(WRF,file="WRF_mslp.RData")
for(i in 1:7) esb.corr(smon[i,],seasons[i],WRF,'MSLP','WRF',ESBW,lat,lon)
