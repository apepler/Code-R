##Comparing all the reanalyses for a single event on June 19 2007

setwd('~/Documents/Data/JE_WRF')
library("RNetCDF")
library("R.matlab")
readMat('~/Documents/Data/Useful.mat')->Useful
mask<-t(Useful$mask)
mask[is.na(mask)]=0
library("akima")

##WRF
f2=open.nc('/srv/ccrc/data23/z3444417/studies/Data/WRF50/2D/WRF50_Lon.nc')
lonW=var.get.nc(f2,'Lon')
f2=open.nc('/srv/ccrc/data23/z3444417/studies/Data/WRF50/2D/WRF50_Lat.nc')
latW=var.get.nc(f2,'Lat')
latt=as.vector(latW)
lont=as.vector(lonW)
lat2=seq(-50,0,0.5)
lon2=seq(100,180,0.5)

fname=paste('/srv/ccrc/data23/z3444417/studies/Data/WRF50/nnrp/R2/psl/raw/monthly/WRF50_R2_slp_2007-06.nc',sep="")
f1=open.nc(fname)
time=var.get.nc(f1,'Times')
date=paste(substr(time,1,4),substr(time,6,7),substr(time,9,10),substr(time,12,13),sep="")
I=which(date==2007061900)
WRF=var.get.nc(f1,'SLP',c(1,1,I),c(215,144,1))
a=as.vector(WRF)
a[which(is.na(a))]=0
b=interp(lont,latt,a,lon2,lat2)
WRF2=b$z
contour(lon2,lat2,WRF2)
par(xpd = NA)
contour(Useful$x,Useful$y,mask,add=T,drawlabels=F)

hh=as.integer(date) %% 100
I=which(hh%%6==0)
date2=date[I]
I=which(date2==2007061900)
f1=open.nc('/srv/ccrc/data23/z3478332/WRF/WRF_slp_2007-06_regrid.nc')
WRFa=var.get.nc(f1,'slp2',c(1,1,I),c(33,19,1))
lata=var.get.nc(f1,'lat2')
lona=var.get.nc(f1,'lon2')
contour(lona,lata,WRFa)
par(xpd = NA)
contour(Useful$x,Useful$y,mask,add=T,drawlabels=F)

WRFb=var.get.nc(f1,'slp1',c(1,1,I),c(54,31,1))
latb=var.get.nc(f1,'lat1')
lonb=var.get.nc(f1,'lon1')
contour(lonb,latb,WRFb)
par(xpd = NA)
contour(Useful$x,Useful$y,mask,add=T,drawlabels=F)

### Okay, ERAI options

f1=open.nc('/srv/ccrc/data23/z3444417/studies/Data/ERAI/global/psl/raw/ERAI_mslp_2000-01_2010-12.nc')
f2=open.nc('/srv/ccrc/data23/z3478332/ERAI/ERAI_mslp_2000-01_2010-12_half.nc')
f3=open.nc('/srv/ccrc/data23/z3478332/ERAI/ERAI_mslp_2000-01_2010-12_regrid.nc')
lat=var.get.nc(f1,'latitude')
lon=var.get.nc(f1,'longitude')
time1=var.get.nc(f1,'time')
time1a=as.Date((time1/24-2),origin="1900-01-03")
I=which(time1a=="2007-06-19")
era=var.get.nc(f1,'msl',c(1,1,I),c(240,120,1),unpack=T)/100
I=which(lat>=(-50) & lat<=0)
lat=lat[I]
J=which(lon>=100 & lon<=180)
lon=lon[J]
era=era[J,I]
contour(lon,lat[33:1],era[,33:1])
par(xpd = NA)
contour(Useful$x,Useful$y,mask,add=T,drawlabels=F)

lat=var.get.nc(f2,'lat')
lon=var.get.nc(f2,'lon')
I=which(time1a=="2007-06-19")
era2=var.get.nc(f2,'msl',c(1,1,I),c(length(lon),length(lat),1))/100
I=which(lat>=(-50) & lat<=0)
lat=lat[I]
J=which(lon>=100 & lon<=180)
lon=lon[J]
era2=era2[J,I]
contour(lon,lat[17:1],era2[,17:1])
par(xpd = NA)
contour(Useful$x,Useful$y,mask,add=T,drawlabels=F)

lat=var.get.nc(f3,'lat')
lon=var.get.nc(f3,'lon')
I=which(time1a=="2007-06-19")
era3=var.get.nc(f3,'msl',c(1,1,I),c(length(lon),length(lat),1))/100
I=which(lat>=(-50) & lat<=0)
lat=lat[I]
J=which(lon>=100 & lon<=180)
lon=lon[J]
era3=era3[J,I]
contour(lon,lat[21:1],era3[,21:1])
par(xpd = NA)
contour(Useful$x,Useful$y,mask,add=T,drawlabels=F)


f1=open.nc('/srv/ccrc/data23/z3478332/WRF/WRF_slp_1979-12_regrid.nc')
lata=var.get.nc(f1,'lat2')
lona=var.get.nc(f1,'lon2')
WRFa=var.get.nc(f1,'slp2')
contour(lona,lata,WRFa)
par(xpd = NA)
contour(Useful$x,Useful$y,mask,add=T,drawlabels=F)

WRFb=var.get.nc(f1,'slp1',c(1,1,I),c(54,31,1))
latb=var.get.nc(f1,'lat1')
lonb=var.get.nc(f1,'lon1')
contour(lonb,latb,WRFb)
par(xpd = NA)
contour(Useful$x,Useful$y,mask,add=T,drawlabels=F)
