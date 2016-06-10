##All years of MSLP averaged across selected months
##Obviously doesn't do cross-year seasons. Yet
rm(list=ls())
setwd('~/Documents/Data/JE_WRF')
library("RNetCDF")
library("abind")
library("akima")
library("R.matlab")
readMat('~/Documents/Data/Useful.mat')->Useful
mask<-t(Useful$mask)
mask[is.na(mask)]=0
source('~/Documents/R/filled.contour3.R')
source('~/Documents/R/filled.legend.R')
source('~/Documents/R/color.palette.R')
pal <- color.palette(c("red","yellow","white","cyan","blue"), c(20,5,5,20))
cm=pal(12)
years=seq(1979,2009)
months=seq(1,12)

NCEP<-array(0,dim=c(144,73,length(years),length(months)))
for(i in 1:length(years))
{
  fname=paste('/srv/ccrc/data23/z3478332/NCEP/slp.',years[i],'.nc',sep="")
  f1=open.nc(fname)
  nt=dim.inq.nc(f1,'time')
  nt<-nt$length
  time=var.get.nc(f1,'time')
  lat2=var.get.nc(f1,'lat')
  lon2=var.get.nc(f1,'lon')
  time2=as.Date((time/24-2),origin="0001-01-01")
  mm=as.numeric(format(time2,"%m"))
  hh=(time %% 24)
  slp=var.get.nc(f1,'slp',c(1,1,1),c(144,73,nt),unpack=T)/100
  for(m in 1:12) NCEP[,,i,m]=apply(slp[,,(hh==0 & mm==m)],c(1,2),mean)
  close.nc(f1)
}
latn=lat2[lat2>(-60) & lat2<10]
lonn=lon2[lon2>100 & lon2<180]
NCEP=NCEP[(lon2>100 & lon2<180),(lat2>(-60) & lat2<10),,]

f2=open.nc('/srv/ccrc/data23/z3444417/Data/WRF/2D/WRF50_Lon.nc')
lon=var.get.nc(f2,'Lon')
f2=open.nc('/srv/ccrc/data23/z3444417/Data/WRF/2D/WRF50_Lat.nc')
lat=var.get.nc(f2,'Lat')

load("~/Documents/Data/JE_WRF/WRF_MSLP.RData")
# WRF<-array(0,dim=c(215,144,length(years),length(months)))
# for(i in 1:length(years))
#   for(m in 1:12)
#   {
#     fname=paste('/srv/ccrc/data23/z3444417/Data/WRF/nnrp/R2/psl/raw/d01/WRF_mslp_R2_d01_',years[i],'-',sprintf("%02d",m),'.nc',sep="")
#     f1=open.nc(fname)
#     nt=dim.inq.nc(f1,'Time')
#     nt<-nt$length
#     slp=var.get.nc(f1,'SLP',c(1,1,1),c(215,144,nt))
#     time=var.get.nc(f1,'Times')
#     time=as.integer(substr(time,12,13))
#     close.nc(f1)
#     WRF[,,i,m]=apply(slp[,,time==0],c(1,2),mean,na.rm=T)
#   }
latt=as.vector(lat)
lont=as.vector(lon)
save(lat,lon,WRF,file="~/Documents/Data/JE_WRF/WRF_MSLP.RData")

NCEP2=apply(NCEP,c(1,2,3),mean,na.rm=T)
WRF2<-array(0,dim=c(length(lonn),length(latn),length(years)))
for(i in 1:length(years))
{
  a=as.vector(apply(WRF[,,i,],c(1,2),mean,na.rm=T))
  b=interp(lont,latt,a,lonn,latn)
  WRF2[,,i]=b$z
}

diff=WRF2-NCEP2
for(i in 1:31)
{
    tiff(file=paste('Plots/slp/MSLP_diff_',years[i],'.tiff',sep=""), height=450, width=600)
    par(plt = c(0.1,0.8,0.1,0.9),las = 1,cex.axis = 1)
    filled.contour3(lonn,latn[27:1],diff[,27:1,i],levels=seq(-3,3,0.5),col=cm)
    par(xpd = NA)
    contour(Useful$x,Useful$y,mask,add=T,drawlabels=F)
    par(new = "TRUE",plt = c(0.85,0.9,0.1,0.9),las = 1,cex.axis = 1)
    filled.legend(lonn,latn[27:1],levels=seq(-3,3,0.5),xlim = c(min(xintercepts),max(xintercepts)),ylim = c(min(slopes),max(slopes)),col=cm,xlab = "",ylab = "")
    dev.off()
}

corN<-corW<-corD<-matrix(0,31,27) ##Corr with year
for(i in 1:31)
  for(j in 1:27)
  {
    corN[i,j]=cor(NCEP2[i,j,],years)
    corW[i,j]=cor(WRF2[i,j,],years)
    corD[i,j]=cor(diff[i,j,],years)
  }

tN<-tW<-tD<-matrix(NaN,31,27)
for(i in 1:31)
  for(j in 1:25)
  {
    a = lm(NCEP2[i,j,]~years)
    tN[i,j] = a$coefficients[[2]]
    a = lm(WRF2[i,j,]~years)
    tW[i,j] = a$coefficients[[2]]
    a = lm(diff[i,j,]~years)
    tD[i,j] = a$coefficients[[2]]
  }

tiff(file='Plots/slp/MSLPtrend_NCEP_ann.tiff', height=600, width=600)
par(plt = c(0.1,0.8,0.1,0.9),las = 1,cex.axis = 1)
filled.contour3(lonn,latn[27:1],tN[,27:1]*10,levels=seq(-0.6,0.6,0.1),col=cm)
par(xpd = NA)
contour(Useful$x,Useful$y,mask,add=T,drawlabels=F)
par(new = "TRUE",plt = c(0.85,0.9,0.1,0.9),las = 1,cex.axis = 1)
filled.legend(lonn,latn[27:1],levels=seq(-0.6,0.6,0.1),col=cm,xlim = c(min(xintercepts),max(xintercepts)),ylim = c(min(slopes),max(slopes)),xlab = "",ylab = "")
dev.off()
tiff(file='Plots/slp/MSLPtrend_WRF_ann.tiff', height=600, width=600)
par(plt = c(0.1,0.8,0.1,0.9),las = 1,cex.axis = 1)
filled.contour3(lonn,latn[27:1],tW[,27:1]*10,levels=seq(-0.6,0.6,0.1),col=cm)
par(xpd = NA)
contour(Useful$x,Useful$y,mask,add=T,drawlabels=F)
par(new = "TRUE",plt = c(0.85,0.9,0.1,0.9),las = 1,cex.axis = 1)
filled.legend(lonn,latn[27:1],levels=seq(-0.6,0.6,0.1),col=cm,xlim = c(min(xintercepts),max(xintercepts)),ylim = c(min(slopes),max(slopes)),xlab = "",ylab = "")
dev.off()
tiff(file='Plots/slp/MSLPtrend_diff_ann.tiff', height=600, width=600)
par(plt = c(0.1,0.8,0.1,0.9),las = 1,cex.axis = 1)
filled.contour3(lonn,latn[27:1],tD[,27:1]*10,levels=seq(-0.6,0.6,0.1),col=cm)
par(xpd = NA)
contour(Useful$x,Useful$y,mask,add=T,drawlabels=F)
par(new = "TRUE",plt = c(0.85,0.9,0.1,0.9),las = 1,cex.axis = 1)
filled.legend(lonn,latn[27:1],levels=seq(-0.6,0.6,0.1),col=cm,xlim = c(min(xintercepts),max(xintercepts)),ylim = c(min(slopes),max(slopes)),xlab = "",ylab = "")
dev.off()

##Diffs and correlations - by season

# seasons=c("ann","warm","cool","mam","jja","son","djf")
# snum=cbind(c(1,11,5,3,6,9,12),c(12,4,10,5,8,11,2))
# for(s in 1:7)
# {
#   if(snum[s,2]>snum[s,1])
#   {
#     WRF2<-NCEP2<-array(0,dim=c(length(lonn),length(latn),length(years)))
#     for(i in 1:length(years))
#     {
#       a=as.vector(apply(WRF[,,i,snum[s,1]:snum[s,2]],c(1,2),mean,na.rm=T))
#       b=interp(lont,latt,a,lonn,latn)
#       WRF2[,,i]=b$z
#       NCEP2[,,i]=apply(NCEP[,,i,snum[s,1]:snum[s,2]],c(1,2),mean,na.rm=T)
#     }
#   }
#   else
#   {
#     WRF2<-NCEP2<-array(0,dim=c(length(lonn),length(latn),(length(years)-1)))
#     for(i in 1:(length(years)-1))
#     {
#       ww=abind(WRF[,,i,snum[s,1]:12],WRF[,,i+1,1:snum[s,2]],along=3)
#       a=as.vector(apply(ww,c(1,2),mean,na.rm=T))
#       b=interp(lont,latt,a,lonn,latn)
#       WRF2[,,i]=b$z
#       nn=abind(NCEP[,,i,snum[s,1]:12],NCEP[,,i+1,1:snum[s,2]],along=3)
#       NCEP2[,,i]=apply(nn,c(1,2),mean,na.rm=T)
#     }
#   }
#   
#   fout1=paste('Plots/slp/MSLP_NCEP_',seasons[s],'.tiff',sep="")
#   fout2=paste('Plots/slp/MSLP_WRF_',seasons[s],'.tiff',sep="")  
#   fout3=paste('Plots/slp/MSLP_diff_',seasons[s],'.tiff',sep="")
#   tiff(file=fout1, height=600, width=600)
#   contour(lonn,latn[27:1],apply(NCEP2[,27:1,],c(1,2),mean),levels=seq(970,1030,2))
#   contour(Useful$x,Useful$y,mask,add=T,drawlabels=F)
#   dev.off()
#   tiff(file=fout2, height=600, width=600)
#   contour(lonn,latn[27:1],apply(WRF2[,27:1,],c(1,2),mean),levels=seq(970,1030,2))
#   contour(Useful$x,Useful$y,mask,add=T,drawlabels=F)
#   dev.off()
#   ##And now, plot the difference!
#   tiff(file=fout3, height=450, width=600)
#   par(plt = c(0.1,0.8,0.1,0.9),las = 1,cex.axis = 1)
#   filled.contour3(lonn,latn[27:1],apply(WRF2[,27:1,],c(1,2),mean)-apply(NCEP2[,27:1,],c(1,2),mean),levels=seq(-3,3,0.5),col=cm)
#   par(xpd = NA)
#   contour(Useful$x,Useful$y,mask,add=T,drawlabels=F)
#   par(new = "TRUE",plt = c(0.85,0.9,0.1,0.9),las = 1,cex.axis = 1)
#   filled.legend(lonn,latn[27:1],levels=seq(-3,3,0.5),xlim = c(min(xintercepts),max(xintercepts)),ylim = c(min(slopes),max(slopes)),col=cm,xlab = "",ylab = "")
#   dev.off()
# }

##Let's see if we can repeat somewhat for upper level data
##Start with WRF...

f1=open.nc('/srv/ccrc/data23/z3478332/WRF/WRF_gph.nc')
latW=var.get.nc(f1,'lat')
lonW=var.get.nc(f1,'lon')
WRF=var.get.nc(f1,'Z')
time=var.get.nc(f1,"time")
yy=floor(time/100)
mm=time%%100

f2=open.nc('/srv/ccrc/data23/z3478332/NCEP/hgt.mon.mean.nc')
latN=var.get.nc(f2,'lat')
lonN=var.get.nc(f2,'lon')
NCEP=var.get.nc(f2,'hgt',c(1,1,2,373),c(144,73,5,372)) ##Only 925/850/700/500
NCEP=NCEP[,,-4,]

latn=latN[latN>(-60) & latN<10]
lonn=lonN[lonN>100 & lonN<180]
NCEP=NCEP[(lonN>100 & lonN<180),(latN>(-60) & latN<10),,]
latt=as.vector(latW)
lont=as.vector(lonW)

seasons=c("ann","warm","cool","mam","jja","son","djf")
snum=cbind(c(1,11,5,3,6,9,12),c(12,4,10,5,8,11,2))
levels=c(925,850,700,500)
for(s in 1:7)
  for(h in 1:4)
  {
    if(snum[s,2]>snum[s,1])
    {
      WRF2<-NCEP2<-array(0,dim=c(31,27,length(years)))
      for(i in 1:length(years))
      {
        I=which(yy==years[i] & mm >=snum[s,1] & mm<=snum[s,2])
        a=as.vector(apply(WRF[,,h,I],c(1,2),mean,na.rm=T))
        J=which(!is.na(a))
        b=interp(lont[J],latt[J],a[J],lonn,latn)
        WRF2[,,i]=b$z
        NCEP2[,,i]=apply(NCEP[,,h,I],c(1,2),mean,na.rm=T)
      }
    }
    else
    {
      WRF2<-NCEP2<-array(0,dim=c(31,27,(length(years)-1)))
      for(i in 1:(length(years)-1))
      {
        I=which((yy==years[i] & mm >=snum[s,1]) | (yy==years[i+1] & mm<=snum[s,2]))
        a=as.vector(apply(WRF[,,h,I],c(1,2),mean,na.rm=T))
        J=which(!is.na(a))
        b=interp(lont[J],latt[J],a[J],lonn,latn)
        WRF2[,,i]=b$z
        NCEP2[,,i]=apply(NCEP[,,h,I],c(1,2),mean,na.rm=T)
      }
    }
    fout1=paste('Plots/gph/GPH_',levels[h],'_NCEP_',seasons[s],'.tiff',sep="")
    fout2=paste('Plots/gph/GPH_',levels[h],'_WRF_',seasons[s],'.tiff',sep="")  
    fout3=paste('Plots/gph/GPH_',levels[h],'_diff_',seasons[s],'.tiff',sep="")
    tiff(file=fout1, height=600, width=600)
    contour(lonn,latn[27:1],apply(NCEP2[,27:1,],c(1,2),mean),levels=seq(0,8000,20))
    contour(Useful$x,Useful$y,mask,add=T,drawlabels=F)
    dev.off()
    tiff(file=fout2, height=600, width=600)
    contour(lonn,latn[27:1],apply(WRF2[,27:1,],c(1,2),mean),levels=seq(0,8000,20))
    contour(Useful$x,Useful$y,mask,add=T,drawlabels=F)
    dev.off()
    ##And now, plot the difference!
    tiff(file=fout3, height=450, width=600)
    par(plt = c(0.1,0.8,0.1,0.9),las = 1,cex.axis = 1)
    filled.contour3(lonn,latn[27:1],apply(WRF2[,27:1,],c(1,2),mean)-apply(NCEP2[,27:1,],c(1,2),mean),levels=seq(-30,30,5),col=cm)
    par(xpd = NA)
    contour(Useful$x,Useful$y,mask,add=T,drawlabels=F)
    par(new = "TRUE",plt = c(0.85,0.9,0.1,0.9),las = 1,cex.axis = 1)
    filled.legend(lonn,latn[27:1],levels=seq(-30,30,5),xlim = c(min(xintercepts),max(xintercepts)),ylim = c(min(slopes),max(slopes)),col=cm,xlab = "",ylab = "")
    dev.off()  
  }

